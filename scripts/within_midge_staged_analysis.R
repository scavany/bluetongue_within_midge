## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspace
rm(list = ls())

## install necessary packages
if(!require(deSolve)){install.packages('deSolve'); library(deSolve)}
if(!require(BayesianTools)){install.packages('BayesianTools'); library(BayesianTools)}

## load within-vector model of BTV infection
source("./within_midge_staged_fn.R")

## Mellor "data"
mellor <- read.csv("../data/Mellor_points.csv",header=F)
names(mellor) <- c("days","titre")
mellor$days <- floor(mellor$days+0.5)
initial.titre <- 10^mellor$titre[1]

## fixed parameters
k <- sqrt(1e7)
epsilon <- 2
eip <- 12
tau <- eip - epsilon
n <- 3

## run model
times = seq(from = 0, to = 130, by = 0.1)
beta <- 0.924
p.m <- 5.74
p.s <- p.m*889 ## p.s > p.m
c <- (p.s + p.m) / k
parameters = c(c.m = c, c.h = c, beta.m = beta, beta.s = beta, n = n,
               p.m = p.m, p.s = p.s, epsilon = epsilon, tau = tau)
V.m <- 1; T.m <- 1; T.s <- 1 # equations are proportions
state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s, rep(0,parameters["n"]+1))
state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                 "V.h","T.s",paste0("I.s",0:(parameters["n"]-1)),"I.sn")
names(state) <- state.names
                 

dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
## plot(dat[,1],1 * (dat[,"V.m"] + dat[,"V.h"]))

## generate plot
tiff("../figures/model_output_fixed_EIP.tif",res=600,
     compression="lzw",height=600*6,width=600*6)
par(mar=c(5.1,4.1,4.1,2.1))
layout(mat = matrix(c(1:8,rep(9,4)), nrow = 3,byrow=TRUE))
plot(dat[,1], dat[,"V.m"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'V.m', las = 1)
plot(dat[,1], dat[,"T.m"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'T.m', las = 1)
plot(dat[,1], dat[,"I.m0"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'I.m0', las = 1)
plot(dat[,1], dat[,"I.mn"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'I.mn', las = 1)
plot(dat[,1], dat[,"V.h"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'V.h', las = 1)
plot(dat[,1], dat[,"T.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'T.s', las = 1)
plot(dat[,1], dat[,"I.s0"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'I.s0', las = 1)
plot(dat[,1], dat[,"I.sn"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'I.sn', las = 1)
plot(dat[,1],initial.titre*(dat[,"V.m"] + dat[,"V.h"]),lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^mellor$titre, initial.titre*(dat[,"V.m"] + dat[,"V.h"]))))
points(mellor$day,10^mellor$titre)
dev.off()

### MCMC
## range of free parameters
p.s.mult <- log(c(1e1,1e5))
p.m <- log(c(1e-3,1e3))
## epsilon <- log(c(1,3))
## eip <- log(c(10,14))
beta <- log(c(1e-3,1e3))
## n.dbl <- log(c(1,11))
## pars <- matrix(c(p.s.mult,p.m,epsilon,eip,beta,n.dbl),nrow=6,byrow=TRUE,
pars <- matrix(c(p.s.mult,p.m,beta),nrow=3,byrow=TRUE,
               dimnames=list(c("p.s.mult","p.m","beta"),
                             c("lower","upper")))

## log-likelihood function
LL <- function(par) {
    p.s.mult <- as.numeric(exp(par[1]))
    p.m <- as.numeric(exp(par[2]))
    ## epsilon <- as.numeric(exp(par[3]))
    ## eip <- as.numeric(exp(par[4]))
    beta <- as.numeric(exp(par[3]))
    ## n.dbl <- floor(as.numeric(exp(par[6])))

    n <- 2
    k <- sqrt(1e7)
    epsilon <- 2
    eip <- 12
    tau <- eip - epsilon
    p.s <- p.s.mult*p.m
    times = seq(from = 0, to = 13, by = 0.1)
    c <- (p.s + p.m) / k
    
    parameters = c(c.m = c, c.h = c, beta.m = beta, beta.s = beta, n = n,
                   p.m = p.m, p.s = p.s, epsilon = epsilon, tau = tau)
    V.m <- 1; T.m <- 1; T.s <- 1
    state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s, rep(0,parameters["n"]+1))
    state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                     "V.h","T.s",paste0("I.s",0:(parameters["n"]-1)),"I.sn")
    names(state) <- state.names
                
    dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
    data.days <- which(dat[,"time"] %in% mellor$days)
    titre.model <- initial.titre*(dat[,"V.m"] + dat[,"V.h"])[data.days]
    titre.data <- floor(10^mellor$titre+0.5)

    return(sum(dpois(titre.data,titre.model,log=TRUE)))
}

## sample posterior using MCMC
bayesianSetup = createBayesianSetup(LL,lower=pars[,"lower"],upper=pars[,"upper"])
settings = list(message = T)
out = runMCMC(bayesianSetup, settings = settings)
samples = getSample(out)
save(out,file="mcmc_out.RData")
exp(colMeans(samples))
plot(out)

###  Fit with optim
## fixed parameters
times = seq(from = 0, to = 13, by = 0.1)
pars.init <- c(n=log(3),beta=log(1),p.m=log(2),p.s.mult=log(1000))

optim.fn <- function(pars, times) {
    k <- sqrt(1e7)
    epsilon <- 2
    eip <- 12
    tau <- eip - epsilon
    n <- floor(as.numeric(exp(pars[1]))+0.5)
    beta <- as.numeric(exp(pars[2]))
    p.m <- as.numeric(exp(pars[3]))
    p.s <- p.m*as.numeric(exp(pars[4])) ## p.s > p.m
    c <- (p.s + p.m) / k
    parameters = c(c.m = c, c.h = c, beta.m = beta, beta.s = beta, n = n,
                   p.m = p.m, p.s = p.s, epsilon = epsilon, tau = tau)
    V.m <- 1; T.m <- 1; T.s <- 1 # equations are proportions
    state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s, rep(0,parameters["n"]+1))
    state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                     "V.h","T.s",paste0("I.s",0:(parameters["n"]-1)),"I.sn")
    names(state) <- state.names

    dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
    data.days <- which(dat[,"time"] %in% mellor$days)
    logtitre.model <- log(initial.titre*(dat[,"V.m"] + dat[,"V.h"])[data.days],10)
    logtitre.data <- mellor$titre
    return(sum(abs(logtitre.model - logtitre.data)))
}

optim.out <- optim(pars.init,optim.fn,times=times)
save(optim.out,file="optim_out.RData")
