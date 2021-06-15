## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspacer
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(BayesianTools,psych,bbmle,deSolve,grDevices)

## load within-vector model of BTV infection
source("./within_midge_staged_fn.R")

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
initial.titre <- 10^fu$log10.mean.titre[1] 

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)

## Set up timing vectors
times = seq(from = 0, to = 14*24, by = 1)
data.hours.it <- which(times %in% fu.intrathoracic.hours)
model.hours.it <- which(fu.intrathoracic.hours %in% times)
titre.data.it <- floor(fu.intrathoracic$titre[model.hours.it]+0.5)
data.hours <- which(times %in% fu$Time.pi..hours.)
model.hours <- which(fu$Time.pi..hours. %in% times)
titre.data <- floor(10^fu$log10.mean.titre[model.hours] + 0.5)

## Fit logistic curve to the proportion positive
prop.positive.fn <- function(times,prop.disseminated,k,t0) {
    1 - (1 - prop.disseminated)/(1 + exp(-k*(times - t0)))
}
logistic.optim.fn <- function(par) {
    prop.positive <- prop.positive.fn(times,par[1],par[2],par[3])
    return(sum((prop.positive[data.hours] - fu$Detection.rate[model.hours])^2))
}
logistic.optim.out <- optim(c(0.4,1,0),logistic.optim.fn)
logistic.optim.out
plot(times,prop.positive.fn(times,logistic.optim.out$par[1],
                            logistic.optim.out$par[2],
                            logistic.optim.out$par[3]),type='l',
     ylim=c(0,1),ylab="Proportion positive")
points(fu$Time.pi..hours.,fu$Detection.rate)
prop.positive.vec <- prop.positive.fn(times,logistic.optim.out$par[1],
                                      logistic.optim.out$par[2],
                                      logistic.optim.out$par[3])
prop.disseminated <- min(prop.positive.vec)

###  Fit with mle2 - do this jointly with secondary stage???
## fixed parameter
n <- 1 # number of incubation steps in midgut

## The NLL function
NLL <- function(pars) {
    ## Secondary tissues, intrathoracic first
    with(as.list(pars),{
        beta.s <- as.numeric(exp(lbeta.s))
        p.s <- as.numeric(exp(lp.s))
        c.s <- p.s / geometric.mean(tail(fu.intrathoracic$titre))
        parameters = c(c.s = c.s, beta.s = beta.s, p.s = p.s)
        V.h <- fu.intrathoracic$titre[1]; T.s <- 1 
        state = c(V.h, T.s, 0)
        state.names <- c("V.h","T.s","I.s")
        names(state) <- state.names
        dat = ode(y = state, times = times, func = wv.BTV.intrathoracic, parms = parameters)
        titre.model.it <- dat[data.hours.it,"V.h"]
        ## Then full infection, accounting for both succesful and failed infections
        beta.m <- as.numeric(exp(lbeta.m))
        p.m <- p.s/1e4#as.numeric(exp(lp.m))
        c.m <- as.numeric(exp(lc.m))
        epsilon <- as.numeric(exp(lepsilon))
        k <- geometric.mean(tail(10^fu$log10.mean.titre))/geometric.mean(tail(fu.intrathoracic$titre)) -
            p.m*(1-exp(-beta.m*initial.titre))/p.s
        ## print(k)
        if (k<0) {return(Inf)}
        parameters = c(c.m = c.m, c.s = c.s, beta.m = beta.m, beta.s = beta.s, n = n,
                       p.m = p.m, p.s = p.s, epsilon = epsilon, k=as.numeric(k))
        V.m <- initial.titre; T.m <- 1; T.s <- 1 # equations are proportions
        state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s,0, 0)
        state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                         "V.h","T.s","I.s0","I.s")
        names(state) <- state.names
        dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
        titre.disseminated <- rowSums(dat[data.hours,c("V.m","V.h")])
        titre.failure <- (V.m*exp(-c.m*times))[data.hours]
        prop.failure <- prop.positive.vec[data.hours] - prop.disseminated
        titre.total <- (titre.disseminated*prop.disseminated + titre.failure * prop.failure) /
            prop.positive.vec[data.hours]
        ## return(-sum(dpois(titre.data,titre.disseminated,log=TRUE)))
        return(-sum(dpois(titre.data,titre.total,log=TRUE))
               -sum(dpois(titre.data.it,titre.model.it,log=TRUE)))
    })
}
LL <- function(par) {
    pars <- log(par)
    names(pars) <- names(pars.init)
    -NLL(pars)
}

pars.init <- c(lbeta.m=log(1e-5),lp.m=log(1e-2),lc.m=log(1e-2),lepsilon=log(2),
               lbeta.s=log(1e-3),lp.s=log(1e2))
optim.out <- optim(pars.init,NLL)
save(optim.out,file="optim_out.RData")
load("optim_out.RData")

lower <- c(lbeta.m=1e-7,lp.m=1e-3,lc.m=1e-4,lepsilon=1)#,
           ##lbeta.s=1e-7,lp.s=1e-3)
upper <- c(lbeta.m=1e-1,lp.m=1e3,lc.m=1e0,lepsilon=3)#,
           ##lbeta.s=1e-1,lp.s=1e3)
bayesianSetup = createBayesianSetup(LL,lower=lower,upper=upper)
mcmc.out = runMCMC(bayesianSetup)#,
                   ##settings=list(iterations=9e5))
    

### Now run the intrathoracic model
beta.s <- as.numeric(exp(optim.out$par["lbeta.s"]))
p.s <- as.numeric(exp(optim.out$par["lp.s"]))
c.s <- as.numeric(p.s / geometric.mean(tail(fu.intrathoracic$titre)))
parameters = c(c.s = c.s, beta.s = beta.s, 
               p.s = p.s)
inoculum.factor <- 1
V.h <- inoculum.factor*fu.intrathoracic$titre[1]; T.s <- 1 
state <- c(V.h,T.s,0) 
names(state) <- c("V.h","T.s","I.s")
dat.intrathoracic = ode(y = state, times = times, func = wv.BTV.intrathoracic,
                        parms = parameters)

## generate plot
tiff(paste0("../figures/model_output_fixed_EIP_Fu_intrathoracic_inoculum_x",inoculum.factor,".tif"),
     res=600,
     compression="lzw",height=600*4,width=600*6)
par(mar=c(5.1,4.1,4.1,2.1))
layout(mat = matrix(c(1:3,rep(4,3)), nrow = 2,byrow=TRUE))
plot(dat.intrathoracic[,1], dat.intrathoracic[,"V.h"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'V.h', las = 1)
plot(dat.intrathoracic[,1], dat.intrathoracic[,"T.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'T.s', las = 1)
plot(dat.intrathoracic[,1], dat.intrathoracic[,"I.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'I.s', las = 1)
plot(dat.intrathoracic[,1],dat.intrathoracic[,"V.h"],lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(fu.intrathoracic$titre, dat.intrathoracic[,"V.h"])))
points(fu.intrathoracic.hours,fu.intrathoracic$titre)
dev.off()

### Now do the two stage model
epsilon <- as.numeric(exp(optim.out$par["lepsilon"]))
beta.m <- as.numeric(exp(optim.out$par["lbeta.m"]))
p.m <- as.numeric(exp(optim.out$par["lp.m"]))
c.m <- as.numeric(exp(optim.out$par["lc.m"]))
k <- geometric.mean(tail(10^fu$log10.mean.titre))/geometric.mean(tail(fu.intrathoracic$titre)) -
    p.m*(1-exp(-beta.m*initial.titre))/p.s
inoculum.factor <- 1
parameters = c(c.m = c.m, c.s = c.s, beta.m = beta.m, beta.s = beta.s, n = n,
               p.m = p.m, p.s = p.s, epsilon = epsilon, k=as.numeric(k))
V.m <- inoculum.factor*initial.titre; T.m <- 1; T.s <- 1 # equations are proportions
state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s,0, 0)
state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                 "V.h","T.s","I.s0","I.s")
names(state) <- state.names
dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
titre.disseminated <- rowSums(dat[,c("V.m","V.h")])

## Get model output for failed infections
titre.failure <- (V.m*exp(-c.m*times))
prop.failure <- prop.positive.vec - prop.disseminated
titre.total <- (titre.disseminated*prop.disseminated + titre.failure * prop.failure) / prop.positive.vec

## generate plot
tiff(paste0("../figures/model_output_fixed_EIP_Fu_inoculum_x",inoculum.factor,".tif"),
            res=600, compression="lzw",height=600*6,width=600*6)
par(mar=c(5.1,4.1,4.1,5.1))
layout(mat = matrix(c(1:6,rep(7,3)), nrow = 3,byrow=TRUE))
plot(dat[,1], dat[,"V.m"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'V.m', las = 1)
plot(dat[,1], dat[,"T.m"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Hours)', ylab = 'T.m', las = 1)
plot(dat[,1], dat[,"I.mn"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Hours)', ylab = 'I.mn', las = 1)
plot(dat[,1], dat[,"V.h"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'V.h', las = 1)
plot(dat[,1], dat[,"T.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'T.s', las = 1)
plot(dat[,1], dat[,"I.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'I.s0', las = 1)
plot(dat[,1],titre.total,lwd = 2, bty = 'n', 
     xlab = 'Time (Hours)', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^0.75,10^fu$log10.mean.titre, dat[,"V.m"] + dat[,"V.h"])))
points(fu$Time.pi..hours,10^fu$log10.mean.titre)
lines(dat[,1], titre.failure,col="red",lwd=2)
lines(dat[,1], titre.disseminated,col="green",lwd=2)
abline(h=10^0.75,lty="dashed")
par(new=TRUE)
plot(times,prop.positive.fn(times,logistic.optim.out$par[1],
                                              logistic.optim.out$par[2],
                                              logistic.optim.out$par[3]),
     type='l',
     ylim=c(0,1),lty='dashed',yaxt="n",xaxt='n',xlab="",ylab="",bty="n")
legend(200,0.82,legend=c("Midgut-only infections","Disseminated infections","All infections"),
       col=c("red","green","black"),lwd=2,bty="n")
axis(4,at=seq(0,1,0.5))
mtext("Proportion positive",4,3)
dev.off()


### boneyard
## optim.fn <- function(pars, times) {
##     ## n <- floor(as.numeric(exp(pars[5]))+0.5)
##     beta.m <- as.numeric(exp(pars[1]))
##     p.m <- as.numeric(exp(pars[2]))
##     c.m <- as.numeric(exp(pars[3]))
##     epsilon <- as.numeric(exp(pars[4]))
##     k <- geometric.mean(tail(10^fu$log10.mean.titre))/geometric.mean(tail(fu.intrathoracic$titre)) -
##         p.m*(1-exp(-beta.m*initial.titre))/p.s
##     if (k<0) {return(Inf)}
##     parameters = c(c.m = c.m, c.s = c.s, beta.m = beta.m, beta.s = beta.s, n = n,
##                    p.m = p.m, p.s = p.s, epsilon = epsilon, k=as.numeric(k))
##     V.m <- initial.titre; T.m <- 1; T.s <- 1 # equations are proportions
##     state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s, 0)
##     state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
##                  "V.h","T.s","I.s")
##     names(state) <- state.names
##     dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
##     data.hours <- which(dat[,"time"] %in% fu$Time.pi..hours.)
##     model.hours <- which(fu$Time.pi..hours. %in% dat[,"time"])
##     logtitre.model <- log((dat[,"V.m"] + dat[,"V.h"])[data.hours],10)
##     logtitre.data <- fu$log10.mean.titre[model.hours]
##     return(sum(abs(logtitre.model - logtitre.data)))
## }
## optim.fn <- function(pars, times) {
##     beta.s <- as.numeric(exp(pars[1]))
##     p.s <- as.numeric(exp(pars[2]))
##     c.s <- p.s / geometric.mean(tail(fu.intrathoracic$titre))
##     parameters = c(c.s = c.s, beta.s = beta.s, p.s = p.s)
##     V.h <- fu.intrathoracic$titre[1]; T.s <- 1 # equations are proportions
##     state = c(V.h, T.s, 0)
##     state.names <- c("V.h","T.s","I.s")
##     names(state) <- state.names
##     dat = ode(y = state, times = times, func = wv.BTV.intrathoracic, parms = parameters)
##     data.hours <- which(dat[,"time"] %in% fu.intrathoracic.hours)
##     model.hours <- which(fu.intrathoracic.hours %in% dat[,"time"])
##     logtitre.model <- log((dat[,"V.h"])[data.hours],10)
##     logtitre.data <- log(fu.intrathoracic$titre[model.hours],10)
##     return(sum(abs(logtitre.model - logtitre.data)))
## }
