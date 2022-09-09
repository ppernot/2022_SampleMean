# setup ####
genSample = TRUE # if FALSE, READ samples
genStats  = TRUE # if FALSE, READ stats

nMC   = 1e4
nBoot = 5e3

gPars    = ErrViewLib::setgPars("publish")
cols     = rep(gPars$cols,2)
cols_tr  = rep(gPars$cols_tr,2)
cols_tr2 = rep(gPars$cols_tr2,2)
ltyp     = c(rep(1,length(gPars$cols)),rep(2,length(gPars$cols)))

## Uncertainty Estimators ####
### Function c4 from Huang2020
c4 = function(N)
  sqrt( 2/(N-1) ) * gamma(N/2) / gamma((N-1)/2)

### Phi function from Cox2017
incgam = function(a,x)
  gsl::gamma_inc(a,x)
phi0 = function(N, alpha)
  sqrt( (N-1)/2 * incgam( (N-3)/2, 1/alpha) / incgam( (N-1)/2, 1/alpha) )
# N=2; q=1; print(phi0(N, 2/(N-1)*q^2)) # Check with Table B1 of Cox2017
phi = function(N, alpha, beta){
  sqrt( (N-1)/2 *
          (incgam((N-3)/2,1/alpha) - incgam((N-3)/2,1/beta)) /
          (incgam((N-1)/2,1/alpha) - incgam((N-1)/2,1/beta))
  )
}
### Model functions
u_fr = function(X, v, N,...)  # Frequentist
  sd(X) / sqrt(N)
u_fru = function(X, v, N,...) # Unbiased frequentist
  u_fr(X, v, N) / c4(N)
u_frun = function(X, v, N, k=3,...) # Approx Unbiased for non-normal
  u_fr(X, v, N) * sqrt((N-1)/(N-1.5-(k-3)/4))
u_NIP = function(X, v ,N,...) # Non informative prior
  sqrt((N - 1) / (N - 3)) * u_fr(X, v, N)
u_char = function(X, v, N,...)  # Characteristic uncertainty Cox2022
  0.5 * qt(0.975, df=N-1) * u_fr(X, v, N)
u_MIP = function (X, v, N,...) { # OHagan2021 mildly informative prior
  ds = N + 2
  s  = sd(X)
  vs = (3 * v + (N - 1) * s^2) / (N + 2)
  sqrt(ds / (ds - 2) * vs / N)
}
u_SIP = function (X, v, N,...) { # OHagan2021 strongly informative prior
  ds = N + 7
  s  = sd(X)
  vs = (8 * v + (N - 1) * s^2) / (N + 7)
  sqrt(ds / (ds - 2) * vs / N)
}
u_Cox = function(X, v, N, p=3,...) { # Cox2019 informative prior
  S = (N - 1) * var(X)
  sig2_min = v / p
  sig2_max = v * p
  alpha = 2 * sig2_max / S
  beta  = 2 * sig2_min / S
  phi(N, alpha, beta) * u_fr(X, v, N)
}
u_Cox0 = function(X, v, N, p=3,...) { # Cox2019 informative prior
  S = (N - 1) * var(X)
  sig2_max = v * p
  alpha = 2 * sig2_max / S
  phi0(N, alpha) * u_fr(X, v, N)
}

## Stats ####
rmsu = function(x)
  sqrt(mean(x^2,na.rm = TRUE))
genStatsFun = function(mu,res,uTrue,nBoot = 1000) {
  umean = lapply(res, function(x){
    sel = is.finite(x)
    S  = x[sel]
    v = mean(S)
    uv = sd(S)/length(S)^0.5
    return(c(v,uv))
  })
  umed = lapply(res, function(x){
    sel = is.finite(x)
    S  = x[sel]
    v = median(S)
    bs = bootstrap::bootstrap(S,nBoot,median)
    uv = sd(bs$thetastar)
    return(c(v,uv))
  })
  urmsu = lapply(res, function(x){
    sel = is.finite(x)
    S  = x[sel]
    v = rmsu(S)
    bs = bootstrap::bootstrap(S,nBoot,rmsu)
    uv = sd(bs$thetastar)
    return(c(v,uv))
  })
  urange = lapply(res, function(x){
    sel = is.finite(x)
    S  = x[sel]
    v = IQR(S)
    bs = bootstrap::bootstrap(S,nBoot,IQR)
    uv = sd(bs$thetastar)
    return(c(v,uv))
  })
  upund = lapply(
    res,
    function(x,y) {
      sel = is.finite(x)
      S = x[sel] < y
      v = mean(S)
      bs = bootstrap::bootstrap(S,nBoot,mean)
      uv = sd(bs$thetastar)
      return(c(v,uv))
    },
    y = uTrue)
  up20 = lapply(
    res,
    function(x,y) {
      sel = is.finite(x)
      S = x[sel] > y/1.2 & x[sel] < y*1.2
      v = mean(S)
      bs = bootstrap::bootstrap(S,nBoot,mean)
      uv = sd(bs$thetastar)
      return(c(v,uv))
    },
    y = uTrue)
  return(
    list(
      umean = umean,
      umed  = umed,
      urmsu = urmsu,
      upund = upund,
      up20  = up20,
      urange = urange
    )
  )
}


## Unit variance distributions ####
Norm = function(N)
  rnorm(N)
T3u = function(N,df=3)
  rt(N,df=3) / sqrt(3)
Unifu = function(N)
  runif(N, -sqrt(3), sqrt(3))
Laplu = function(N,df=1)
  normalp::rnormp(N, p = df) / sqrt(df^(2/df)*gamma(3/df)/gamma(1/df))


## Lists of cases
mtab = c('u_fr','u_fru','u_frun','u_NIP','u_char','u_MIP','u_SIP','u_Cox')
ftab = c('Unifu','Norm','Laplu', 'T3u')
ktab = c(    9/5,     3,      6,   Inf) # Kurtosis
v = 1 # True variance

if(genSample) {
  # Generate random samples ####
  set.seed(123)
  N = 4
  resu = resum = list()
  for (k in seq_along(ftab)) {
    resu[[ftab[k]]] = list()
    fun = get(ftab[k])
    resum[[ftab[k]]] = rep(0,nMC)
    for(j in 1:nMC) {
      S = fun(N) # Random sample
      resum[[ftab[k]]][j] = mean(S)
      for (l in seq_along(mtab)) {
        meth = get(mtab[l])
        resu[[ftab[k]]][[mtab[l]]][j] = meth(S,v,N,k=ktab[k])
      }
    }
  }
  resu4  = resu
  resum4 = resum

  N = 40
  resu = resum = list()
  for (k in seq_along(ftab)) {
    resu[[ftab[k]]] = list()
    fun = get(ftab[k])
    resum[[ftab[k]]] = rep(0,nMC)
    for(j in 1:nMC) {
      S = fun(N) # Random sample
      resum[[ftab[k]]][j] = mean(S)
      for (l in seq_along(mtab)) {
        meth = get(mtab[l])
        resu[[ftab[k]]][[mtab[l]]][j] = meth(S,v,N,k=ktab[k])
      }
    }
  }
  resu40  = resu
  resum40 = resum

  save(resu4,resum4,resu40,resum40,
       file ='Tables/samples.Rda')

} else {
  load(file='Tables/samples.Rda')
}
uTrue = 1/sqrt(4)
uTrue40 = 1/sqrt(40)


# Table 1 ####
fmean4 = lapply(resum4,function(x){
  sel = is.finite(x)
  S  = x[sel]
  ErrViewLib::prettyUnc(mean(S),sd(S)/length(S)^0.5,numDig = 1)
})
fsd4   = lapply(resum4,function(x){
  sel = is.finite(x)
  S  = x[sel]
  v = sd(S)
  bs = bootstrap::bootstrap(S,nBoot,sd)
  uv = sd(bs$thetastar)
  ErrViewLib::prettyUnc(v,uv,numDig = 1)
})
fkurt4 = lapply(resum4,function(x) {
  sel = is.finite(x)
  S  = x[sel]
  v = moments::kurtosis(S)
  bs = bootstrap::bootstrap(S,nBoot,moments::kurtosis)
  uv = sd(bs$thetastar)
  ErrViewLib::prettyUnc(v,uv,numDig = 1)
})
fmean40 = lapply(resum40,function(x){
  sel = is.finite(x)
  S  = x[sel]
  ErrViewLib::prettyUnc(mean(S),sd(S)/length(S)^0.5,numDig = 1)
})
fsd40   = lapply(resum40,function(x){
  sel = is.finite(x)
  S  = x[sel]
  v = sd(S)
  bs = bootstrap::bootstrap(S,nBoot,sd)
  uv = sd(bs$thetastar)
  ErrViewLib::prettyUnc(v,uv,numDig = 1)
})
fkurt40 = lapply(resum40,function(x) {
  sel = is.finite(x)
  S  = x[sel]
  v = moments::kurtosis(S)
  bs = bootstrap::bootstrap(S,nBoot,moments::kurtosis)
  uv = sd(bs$thetastar)
  ErrViewLib::prettyUnc(v,uv,numDig = 1)
})

tab = rbind(
  fmean4,fsd4,fkurt4,
  fmean40,fsd40,fkurt40
)
sink(file='Tables/Table1.tex')
knitr::kable(tab, format = 'latex')
sink()

# Stats ####

if(genStats) {
  tabres4 = tabres40 = list()
  for (k in seq_along(ftab)) {
    f   = ftab[k]
    mu  = resum4[[f]]
    res = resu4[[f]]
    tabres4[[f]] = genStatsFun(mu,res,uTrue,nBoot)
    mu  = resum40[[f]]
    res = resu40[[f]]
    tabres40[[f]] = genStatsFun(mu,res,uTrue40,nBoot)
  }
  save(tabres4, tabres40,
       file ='Tables/stats.Rda')
} else {
  load(file='Tables/stats.Rda')
}


# Table 2 ####
f = 'Norm'
tab = rbind(
  apply(as.data.frame(tabres4[[f]]$umean),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres4[[f]]$umed),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres4[[f]]$urmsu),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres4[[f]]$urange),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres4[[f]]$upund),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres4[[f]]$up20),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres40[[f]]$umean),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres40[[f]]$umed),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres40[[f]]$urmsu),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres40[[f]]$urange),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres40[[f]]$upund),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1)),
  apply(as.data.frame(tabres40[[f]]$up20),2,function(x) ErrViewLib::prettyUnc(x[1],x[2],numDig = 1))
  )
sink(file='Tables/Table2.tex')
knitr::kable(tab,format = 'latex')
sink()

# Fig. 1 ####
png(file = 'Figs/Fig_01.png',
    width = 2*gPars$reso, height = gPars$reso)
par(mfrow = c(1,2),
    mar = c(3,3,1,1),
    tcl = gPars$tcl,
    mgp = gPars$mgp,
    pty = gPars$pty,
    lwd = 2*gPars$lwd,
    cex = gPars$cex)

names(cols) = mtab
names(cols_tr2) = mtab
fun = "Norm"
for(num in c(4,40)) {
  resu = get(paste0('resu',num))
  ymax =  max(unlist(lapply(resu[[fun]],function(x) max(density(x)$y))))
  first = TRUE
  for (l in seq_along(mtab)) {
    meth = mtab[l]
    X = resu[[fun]][[meth]]
    if(first) {
      plot(
        density(X),
        col = cols[meth],
        # lwd = 2,
        pch = NA,
        main = '', #paste0('N = ',N,', Dist = ',fun),
        xlab = 'Uncert. on mean',
        xlim = if(num == 4) c(0.0,1.5) else c(0.1,0.22),
        ylim = c(0,ymax),
        xaxs = 'i',
        yaxs = 'i'
      )
      first = FALSE
    } else {
      lines(density(X),col = cols[meth])
    }
  }
  abline(v=1/sqrt(num), lty = 2, col = cols[1])
  grid(lwd=2)
  legend('topright', bty='n', cex=0.75,
         title = paste0('N = ',num),
         legend = c(mtab,'Truth'),
         col = c(cols,1),
         lty = c(rep(1,length(mtab)),2),
         lwd = 2*gPars$lwd)
  box()
}
dev.off()

# Fig. 2 ####
num = 4
props = c('umean','umed','urmsu','urange','upund','up20')
props_pretty = c('Mean','Median','RMS','IQR','P<','P20')
for(p in seq_along(props)) {
  prop = props[p]
  pname = props_pretty[p]

  png(file = paste0('Figs/Fig_02',letters[p],'.png'),
      width = gPars$reso, height = gPars$reso)
  par(mfrow = c(1,1),
      mar = c(4,3,1,1),
      tcl = gPars$tcl,
      mgp = gPars$mgp,
      pty = gPars$pty,
      lwd = 2*gPars$lwd,
      cex = gPars$cex)

  names(cols) = ftab
  names(cols_tr2) = ftab
  names(ltyp) = ftab
  first = TRUE
  for(k in seq_along(ftab)) {
    fun = ftab[k]
    Y = uY = c()
    for(j in seq_along(mtab)) {
      Y[j] = tabres4[[fun]][[prop]][[j]][1]
      uY[j] = tabres4[[fun]][[prop]][[j]][2]
      if(fun == "T3u" & mtab[j] == "u_frun") { # Filter out invalid case
        Y[j]  = NA
        uY[j] = 0
      }
    }

    X = 1:length(Y)
    if(first) {
      plot(
        X,Y,
        type = 'b',
        pch  = 16,
        lty  = ltyp[fun],
        col  = cols[fun],
        main = '',
        xlab = '',
        ylab = pname,
        ylim = if(prop == 'upund' | prop == 'up20') c(0,1)
               else if(prop == 'urange') c(0,0.6)
               else c(0.3,0.9),
        xaxt = 'n',
        xaxs = 'i',
        yaxs = 'i'
      )
      axis(1,at=X,labels = mtab, las=2)
      grid()
      if(!prop %in% c('urange','up20'))
        abline(h = if(prop == 'upund') 0.5 else uTrue,
               lty = 2, col = cols[1])
      first = FALSE
    } else {
      lines(X,Y, type = 'b', pch=16, lty = ltyp[fun], col = cols[fun])
    }
    segments(X,Y-2*uY,X,Y+2*uY,col = cols[fun])

    legend('bottom', ncol = 2, bty = 'n', cex = 0.75,
           # title = paste0('N = ',num),
           legend = ftab,
           col = cols,
           lty = ltyp,
           lwd = 2*gPars$lwd)
    box()
  }
  dev.off()
}

# Fig. 3 ####
num = 40
for(p in seq_along(props)) {
  prop = props[p]
  pname = props_pretty[p]

  png(file = paste0('Figs/Fig_03',letters[p],'.png'),
      width = gPars$reso, height = gPars$reso)
  par(mfrow = c(1,1),
      mar = c(4,3,1,1),
      tcl = gPars$tcl,
      mgp = gPars$mgp,
      pty = gPars$pty,
      lwd = 2*gPars$lwd,
      cex = gPars$cex)

  names(cols) = ftab
  names(cols_tr2) = ftab
  first = TRUE
  for(k in seq_along(ftab)) {
    fun = ftab[k]
    Y = uY = c()
    for(j in seq_along(mtab)) {
      Y[j] = tabres40[[fun]][[prop]][[j]][1]
      uY[j] = tabres40[[fun]][[prop]][[j]][2]
      if(fun == "T3u" & mtab[j] == "u_frun") { # Filter out invalid case
        Y[j]  = NA
        uY[j] = 0
      }
    }
    X = 1:length(Y)
    if(first) {
      plot(
        X,Y,
        type = 'b',
        pch  = 16,
        lty  = ltyp[fun],
        col  = cols[fun],
        main = '',
        xlab = '',
        ylab = pname,
        ylim =if(prop == 'upund' | prop == 'up20') c(0,1)
               else if(prop == 'urange') c(0,0.05)
               else c(0.13,0.17),
        xaxt = 'n',
        xaxs = 'i',
        yaxs = 'i'
      )
      axis(1,at=X,labels = mtab, las=2)
      grid()
      if(!prop %in% c('urange','up20'))
        abline(h = if(prop == 'upund') 0.5 else uTrue40,
               lty = 2, col = cols[1])
      first = FALSE
    } else {
      lines(X,Y, type = 'b', pch=16, lty = ltyp[fun], col = cols[fun])
    }
    segments(X,Y-2*uY,X,Y+2*uY,col = cols[fun])

    legend('bottomleft', ncol = 2, bty = 'n', cex = 0.75,
           # title = paste0('N = ',num),
           legend = ftab,
           col = cols,
           lty = ltyp,
           lwd = 2*gPars$lwd)
    box()
  }
  dev.off()
}

