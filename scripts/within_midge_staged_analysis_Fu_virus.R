## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspacer
rm(list = ls())

## set seed
set.seed(22)

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(BayesianTools,psych,bbmle,deSolve,grDevices,data.table)

## load within-vector model of BTV infection
source("./within_midge_staged_virus_only_fn.R")

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
fu$log10.mean.titre <- log(log(2) * 10^fu$log10.mean.titre,10) # convert to pfu
fu.hours <- fu$Time.pi..hours.
V.m0 <- 10^fu$log10.mean.titre[1]
V.s0 <- 0
V.sinf <- mean(tail(10^fu$log10.mean.titre,10))
titre.data <- floor(10^fu$log10.mean.titre + 0.5)

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic$titre <- fu.intrathoracic$titre * log(2) #convert to pfu
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)
V.s0.it <- fu.intrathoracic$titre[1]
V.sinf.it <- mean(tail(fu.intrathoracic$titre,8))
titre.data.it <- floor(fu.intrathoracic$titre+0.5)

## Try a logistic growth model
times <- seq(1,500,1)
plot(out,log='y',ylim=range(fu.intrathoracic$titre))
parms <- c(K.s=V.sinf.it,g.s=10^-1.3)
state <- c(V.s=V.s0.it)
## out <- ode(state,times,wv.BTV.secondary.logistic,parms)
## plot(out,log='',ylim=range(fu.intrathoracic$titre))
## points(fu.intrathoracic.hours,fu.intrathoracic$titre)
## Looks okay. Maybe use this as the secondary model?

## Fit logistic curve to the proportion positive
prop.positive.fn <- function(times,prop.disseminated,k,t0) {
    1 - (1 - prop.disseminated)/(1 + exp(-k*(times - t0)))
}
logistic.optim.fn <- function(par) {
    prop.positive <- prop.positive.fn(fu.hours,par[1],par[2],par[3])
    return(sum((prop.positive - fu$Detection.rate)^2))
}
logistic.optim.out <- optim(c(0.4,1,0),logistic.optim.fn)
logistic.optim.out
tiff("../figures/logistic_positivity.tif",
     res=600,
     compression="lzw",height=600*6,width=600*6)
plot(seq(0,240,1),prop.positive.fn(seq(0,240,1),logistic.optim.out$par[1],
                            logistic.optim.out$par[2],
                            logistic.optim.out$par[3]),type='l',
     ylim=c(0,1),ylab="Proportion positive",xlab="Time (hours)")
points(fu.hours,fu$Detection.rate)
dev.off()
prop.positive.vec <- prop.positive.fn(fu.hours,logistic.optim.out$par[1],
                                      logistic.optim.out$par[2],
                                      logistic.optim.out$par[3])
prop.disseminated <- min(prop.positive.vec)

###  Fit with optim - do all steps jointly
## The NLL function
NLL <- function(pars) {
    with(as.list(pars),{
        ## Midgut infection
        c.m <- as.numeric(exp(lc.m))
        titre.midgut <- wv.BTV.midgut(fu.hours,c.m,V.m0)
        ## alpha <- -log(1 - prop.disseminated) * c.m / V.m0
        alpha <- as.numeric(exp(lalpha))
        ## Secondary tissues, intrathoracic first
        state <- c(V.m=0,V.s=V.s0.it)
        parms <- c(K.s=V.sinf.it,g.s=g.s,alpha=alpha,c.m=c.m)
        times <- seq(0,500,1)
        titre.model.it <- ode(state,times,wv.BTV.secondary.logistic,
                              parms)[times %in% fu.intrathoracic.hours,"V.s"]
        
        ## Then full infection, accounting for both successful and failed infections
        ## estimate expected time of passing MIB (but try and get an analytic result for this)
        ## times.temp <- seq(0,1e3,0.001)
        ## t.passMIB <- sum(times.temp*d1passMIB(times.temp,
        ##                                  MIBpassrate,V.m0,c.m)*0.001) - 0.5*sum(range(times.temp)*d1passMIB(range(times.temp),MIBpassrate,V.m0,c.m)*0.001)
        ## times.secondary <- times - t.passMIB
        parms <- c(K.s=V.sinf,g.s=g.s,alpha=alpha,c.m=c.m)
        state <- c(V.m=V.m0,V.s=0)
        titre.secondary <- ode(state,times,wv.BTV.secondary.logistic,
                               parms)[times %in% fu.hours,"V.s"]
        titre.disseminated <- titre.midgut + titre.secondary
        prop.failure <- prop.positive.vec - prop.disseminated
        titre.total <- (titre.disseminated * prop.disseminated + titre.midgut * prop.failure) /
            prop.positive.vec
        return(-sum(dpois(titre.data,titre.total,log=TRUE))
               -sum(dpois(titre.data.it,titre.model.it,log=TRUE)))
    })
}
## NLL.lag <- function(pars) {
##     with(as.list(pars),{
##         lag <- as.numeric(exp(llag))
##         ## Secondary tissues, intrathoracic first
##         parms <- c(K.s=V.sinf.it,g.s=g.s)
##         titre.model.it <- ode(state,times,wv.BTV.secondary.logistic,
##                               parms)[times %in% fu.intrathoracic.hours,"V.s"]

##         ## Then full infection, accounting for both successful and failed infections
##         c.m <- as.numeric(exp(lc.m))
##         titre.midgut <- wv.BTV.midgut(fu.hours,c.m,V.m0)
##         MIBpassrate <- -log(1 - prop.disseminated) * c.m / V.m0
##         ## estimate expected time of passing MIB (but try and get an analytic result for this)
##         times.temp <- seq(0,1e3,0.001)
##         t.passMIB <- sum(times.temp*d1passMIB.lag(times.temp,
##                                          MIBpassrate,V.m0,c.m,lag)*0.001) - 0.5*sum(range(times.temp)*d1passMIB.lag(range(times.temp),MIBpassrate,V.m0,c.m,lag)*0.001)
##         times.secondary <- times - t.passMIB
##         parms <- c(K.s=V.sinf,g.s=g.s)
##         titre.secondary <- ode(c(V.s=1),times,wv.BTV.secondary.logistic,
##                                parms)[times.secondary %in% fu.hours,"V.s"]
##         titre.disseminated <- titre.midgut + titre.secondary
##         prop.failure <- prop.positive.vec - prop.disseminated
##         titre.total <- (titre.disseminated * prop.disseminated + titre.midgut * prop.failure) /
##             prop.positive.vec
##         return(-sum(dpois(titre.data,titre.total,log=TRUE))
##                -sum(dpois(titre.data.it,titre.model.it,log=TRUE)))
##     })
## }

pars.init <- c(lc.m=0,g.s=0,lalpha=0)
##pars.init.lag <- c(lc.m=0,g.s=0,llag=log(10))
optim.out <- optim(pars.init,NLL)
##optim.out.lag <- optim(pars.init.lag,NLL.lag)
##2*(length(optim.out$par) - optim.out$value)
##2*(length(optim.out.lag$par) - optim.out.lag$value)
##save(optim.out,optim.out.lag,file="optim_out_virus.RData")
save(optim.out,file="optim_out_virus.RData")

load("optim_out_virus.RData")

### Now plot everything
## Lag term is basically zero so just plot the unlagged model
## lag <- as.numeric(exp(optim.out.lag$par['llag']))
times <- seq(0,max(c(fu.hours,fu.intrathoracic.hours))+24,1)
c.m <- as.numeric(exp(optim.out$par['lc.m']))
g.s <- as.numeric(optim.out$par['g.s'])
alpha <- as.numeric(exp(optim.out$par['lalpha']))#-log(1 - prop.disseminated) * c.m / V.m0
parms <- c(K.s=V.sinf.it,g.s=g.s,alpha=alpha,c.m=c.m)
state <- c(V.m=0,V.s=V.s0.it)
titre.model.it <- ode(state,times,wv.BTV.secondary.logistic,
                      parms)[,"V.s"]
titre.midgut <- wv.BTV.midgut(times,c.m,V.m0)
## ## With lags
## c.s.lag <- as.numeric(exp(optim.out.lag$par['lc.s']))
## titre.model.it.lag <- wv.BTV.secondary(times,c.s.lag,V.s0.it,V.sinf.it)
## c.m.lag <- as.numeric(exp(optim.out.lag$par['lc.m']))
## titre.midgut.lag <- wv.BTV.midgut(times,c.m.lag,V.m0)
## MIBpassrate <- -log(1 - prop.disseminated) * c.m.lag / V.m0
## estimate expected time of passing MIB (but try and get an analytic result for this)
## times.temp <- seq(0,1e3,0.001)
## t.passMIB <- sum(times.temp*d1passMIB(times.temp,
##                                       MIBpassrate,V.m0,c.m)*0.001) - 0.5*sum(range(times.temp)*d1passMIB(range(times.temp),
##                                                                                                          MIBpassrate,V.m0,c.m)*0.001)
## t.passMIB.lag <- sum(times.temp*d1passMIB.lag(times.temp,
##                                          MIBpassrate,V.m0,c.m.lag,lag)*0.001) - 0.5*sum(range(times.temp)*d1passMIB.lag(range(times.temp),
##                                                                                                         MIBpassrate,V.m0,c.m.lag,lag)*0.001)
state <- c(V.m=V.m0,V.s=0)
parms <- c(K.s=V.sinf,g.s=g.s,alpha=alpha,c.m=c.m)
titre.secondary <- ode(state,times,wv.BTV.secondary.logistic,
                       parms)[,"V.s"]
## titre.secondary.lag <- wv.BTV.secondary(times - t.passMIB.lag,c.s.lag,V.s0,V.sinf)
titre.disseminated <- titre.midgut + titre.secondary
## titre.disseminated.lag <- titre.midgut.lag + titre.secondary.lag
prop.positive.vec <-  prop.positive.fn(times,logistic.optim.out$par[1],
                                       logistic.optim.out$par[2],
                                       logistic.optim.out$par[3])
prop.failure <- prop.positive.vec - prop.disseminated
titre.total <- (titre.disseminated * prop.disseminated + titre.midgut * prop.failure) /
    prop.positive.vec
## titre.total.lag <- (titre.disseminated.lag * prop.disseminated + titre.midgut.lag * prop.failure) /
##     prop.positive.vec

## generate plot
tiff("../figures/model_output_virus_only_model.tif",
     res=600,
     compression="lzw",height=600*8,width=600*12)
par(mar=c(5.1,4.1,4.1,5.1),mfrow=c(2,1))
plot(times, titre.model.it, type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (hours)', ylab = 'Viral load (pfu/midge)', las = 1,
     ylim = range(c(titre.model.it,fu.intrathoracic$titre)),log='y')
points(fu.intrathoracic.hours,fu.intrathoracic$titre)
## lines(times,titre.model.it.lag,lty="dashed",lwd=2)
plot(times, titre.total, type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (hours)', ylab = 'Viral load (pfu/midge)', las = 1,
     ylim = range(c(titre.total,10^fu$log10.mean.titre)),log='y')
lines(times,titre.midgut,col='red')
lines(times,titre.secondary,col='green')
## lines(times,titre.total.lag,lty="dashed",lwd=2)
## lines(times,titre.secondary.lag,col='green',lty="dashed")
points(fu.hours,10^fu$log10.mean.titre)
par(new=TRUE)
plot(times,prop.positive.fn(times,logistic.optim.out$par[1],
                            logistic.optim.out$par[2],
                            logistic.optim.out$par[3]),type='l',
     ylim=c(0,1),lty='dashed',yaxt="n",xaxt='n',xlab="",ylab="",bty="n")
axis(4,at=seq(0,1,0.5))
mtext("Proportion positive",4,3)
dev.off()

### Step 1. Draw the wait (calculated dist) and the equilibrium (lnorm) for each midge
## Waits and passes
times <- seq(0,100,0.01)
n <- 1e3
waits <- sample(times[-1],n,replace=TRUE,
                prob=p1passMIB(times[-1],MIBpassrate,V.m0,c.m)-p1passMIB(times[-length(times)],MIBpassrate,V.m0,c.m))
passes <- sample(c(TRUE,FALSE),n,replace=TRUE,prob=c(prop.disseminated,1-prop.disseminated))
## Equilibriums
median <- 10^3.125
meanlog <- log(median)
sdlog <- (5.5-0.75)*log(10)/3.92
equilibriums <- rlnorm(n,meanlog,sdlog)

## Calc traj for each Midge
times <- seq(0,max(c(fu.hours,fu.intrathoracic.hours))+48,1)
titre.midgut <- wv.BTV.midgut(times,c.m,V.m0)
titre.secondary.mat <- matrix(data=NA,ncol=length(passes),nrow=length(times))
for (ii in 1:length(times)) {
    titre.secondary.mat[ii,] <- ifelse(passes,pmax(wv.BTV.secondary(times[ii] - waits,c.s,V.s0,equilibriums),0),0)
}
titre.total <- titre.midgut + titre.secondary.mat

## Plot a subset of the trajectories
tiff("../figures/model_output_example_midge_trajectories.tif",
     res=600,
     compression="lzw",height=600*4,width=600*6)
plot(fu.hours,10^fu$log10.mean.titre,ylim=10^range(c(fu$lower,fu$upper)),log='y',bty="n",
     xlab="Time (hours)",ylab="Viral load (pfu)")
for (jj in sample(n,25)) {
    lines(times,titre.total[,jj])
}
lines(fu.hours,10^fu$upper,lty="dashed")
lines(fu.hours,10^fu$lower,lty="dashed")
dev.off()

## Save output
save(titre.total, titre.midgut,waits, passes, equilibriums,n,file="virus_only_model_basic.RData")
