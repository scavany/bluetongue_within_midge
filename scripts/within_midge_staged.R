# set working directory
setwd('~/Documents/bluetongue_project/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(deSolve)){install.packages('deSolve'); library(deSolve)}

# specify within-vector model of BTV infection
wv.BTV = function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    dV.m = -c.m * V.m
    dT.m = -beta.m * T.m * V.m
    dI.m0 = beta.m * T.m * V.m - n * I.m0 / epsilon
    dI.m <- vector(mode="numeric",length=n-1)
    for (i in 1:(n-1)) {
        diffval <- get(paste0("I.m",i-1)) - get(paste0("I.m",i))
        dI.m[i] <- n * diffval / epsilon
    }
    dI.mn = n * get(paste0("I.m",n-1)) / epsilon
    ## dV.h = p * I.mn - c.h * V.h # only midgut infected cells produce virus
    dV.h = p.m * I.mn + p.s * I.sn - c.h * V.h # all infected cells produce virus
    dT.s = -beta.s * T.s * V.h
    dI.s0 = beta.s * T.s * V.h - n * I.s0 / tau
    dI.s <- vector(mode="numeric",length=n-1)
    for (i in 1:(n-1)) {
        diffval <- get(paste0("I.s",i-1)) - get(paste0("I.s",i))
        dI.s[i] <- n * diffval / tau
    }
    dI.sn = n * get(paste0("I.s",n-1)) / tau
    list(c(dV.m, dT.m, dI.m0, dI.m, dI.mn,
           dV.h, dT.s, dI.s0, dI.s, dI.sn))
  })
}

## fixed parameters
k <- sqrt(1e7)
epsilon <- 2
eip <- 12
tau <- eip - epsilon
n <- 10

## run model
times = seq(from = 0, to = 30, by = 0.01)
beta <- 0.1
p.m <- 10
p.s <- 100*p.m # p.s > p.m
c <- (p.s + p.m) / k
parameters = c(c.m = c, c.h = c, beta.m = beta, beta.s = beta, n = n,
               p.m = p.m, p.s = p.s, epsilon = epsilon, tau = tau)
V.m <- 1; T.m <- 1; T.s <- 1 #equations are proportions
state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s, rep(0,parameters["n"]+1))
state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                 "V.h","T.s",paste0("I.s",0:(parameters["n"]-1)),"I.sn")
names(state) <- state.names
                 

dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
## plot(dat[,1],1 * (dat[,"V.m"] + dat[,"V.h"]))

# generate plot
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
plot(dat[,1],1*(dat[,"V.m"] + dat[,"V.h"]),lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Viral load', las = 1,log="y",type='l')
