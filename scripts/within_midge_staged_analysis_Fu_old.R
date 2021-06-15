## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspacer
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(BayesianTools,psych,bbmle)

## load within-vector model of BTV infection
source("./within_midge_staged_fn.R")

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
initial.titre <- 10^fu$log10.mean.titre[1]

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)

## run model
times = seq(from = 0, to = 10, by = 0.01)
beta.s <- exp(-8.076)
p.s <- exp(10.44)
c.s <- p.s / geometric.mean(tail(fu.intrathoracic$titre))
parameters = c(c.s = c.s, beta.s = beta.s, 
               p.s = p.s)
inoculum.factor <- 1
V.h <- inoculum.factor*fu.intrathoracic$titre[1]; T.s <- 1 # equations are proportions
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
     xlab = 'Time (Days)', ylab = 'V.h', las = 1)
plot(dat.intrathoracic[,1], dat.intrathoracic[,"T.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'T.s', las = 1)
plot(dat.intrathoracic[,1], dat.intrathoracic[,"I.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'I.s', las = 1)
plot(dat.intrathoracic[,1],dat.intrathoracic[,"V.h"],lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(fu.intrathoracic$titre, dat.intrathoracic[,"V.h"])))
points(fu.intrathoracic$day,fu.intrathoracic$titre)
dev.off()

###  Fit with optim
## fixed parameters
times = seq(from = 0, to = 10, by = 0.01)
pars.init <- c(beta=log(1),p.s=log(10))

##fu time points
fu.intrathoracic.days <- floor(fu.intrathoracic$day*10 + 0.5)/10

optim.fn <- function(pars, times) {
    beta.s <- as.numeric(exp(pars[1]))
    p.s <- as.numeric(exp(pars[2]))
    c.s <- p.s / geometric.mean(tail(fu.intrathoracic$titre))
    parameters = c(c.s = c.s, beta.s = beta.s, p.s = p.s)
    V.h <- fu.intrathoracic$titre[1]; T.s <- 1 # equations are proportions
    state = c(V.h, T.s, 0)
    state.names <- c("V.h","T.s","I.s")
    names(state) <- state.names
    dat = ode(y = state, times = times, func = wv.BTV.intrathoracic, parms = parameters)
    data.days <- which(dat[,"time"] %in% fu.intrathoracic.days)
    model.days <- which(fu.intrathoracic.days %in% dat[,"time"])
    logtitre.model <- log((dat[,"V.h"])[data.days],10)
    logtitre.data <- log(fu.intrathoracic$titre[model.days],10)
    return(sum(abs(logtitre.model - logtitre.data)))
}

NLL <- function(pars, times) {
    beta.s <- as.numeric(exp(pars[1]))
    p.s <- as.numeric(exp(pars[2]))
    c.s <- p.s / geometric.mean(tail(fu.intrathoracic$titre))
    parameters = c(c.s = c.s, beta.s = beta.s, p.s = p.s)
    V.h <- fu.intrathoracic$titre[1]; T.s <- 1 # equations are proportions
    state = c(V.h, T.s, 0)
    state.names <- c("V.h","T.s","I.s")
    names(state) <- state.names
    dat = ode(y = state, times = times, func = wv.BTV.intrathoracic, parms = parameters)
    data.days <- which(dat[,"time"] %in% fu.intrathoracic.days)
    model.days <- which(fu.intrathoracic.days %in% dat[,"time"])
    titre.model <- dat[data.days,"V.h"]
    titre.data <- floor(fu.intrathoracic$titre[model.days]+0.5)
    return(-sum(dpois(titre.data,titre.model,log=TRUE)))
}    

optim.out.intrathoracic <- optim(pars.init,NLL,times=times)
save(optim.out.intrathoracic,file="optim_out_intrathoracic.RData")

### Now do the two stage model, with hemocoel parameters as above
## Load optim results
load("./optim_out_intrathoracic.RData",verbose=TRUE)

## fixed parameters
epsilon <- exp(0.3665)
n <- 1
prop.disseminated <- 0.354

## Fit logistic curve to the proportion positive
times = seq(from = 0, to = 14*24, by = 1)
prop.positive.fn <- function(times,prop.disseminated,k,t0) {
    1 - (1 - prop.disseminated)/(1 + exp(-k*(times - t0)))
}
logistic.optim.fn <- function(par) {
    prop.positive <- prop.positive.fn(times,par[1],par[2],par[3])
    data.days <- which(times %in% fu$Time.pi..hours.)
    model.days <- which(fu$Time.pi..hours. %in% times)
    return(sum((prop.positive[data.days] - fu$Detection.rate[model.days])^2))
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
points(times,prop.positive.vec,pch=10)

## run model - disseminated infections
beta.m <- exp(0.1755)
beta.s <- as.numeric(exp(optim.out.intrathoracic$par["beta"])/24)
p.m <- exp(3.403)
p.s <- as.numeric(exp(optim.out.intrathoracic$par["p.s"])/24)
c.m <- exp(-2.082)
c.s <- as.numeric(p.s / geometric.mean(tail(fu.intrathoracic$titre)))
k <- geometric.mean(tail(10^fu$log10.mean.titre))/geometric.mean(tail(fu.intrathoracic$titre)) -
    p.m*(1-exp(-beta.m*initial.titre))/p.s
inoculum.factor <- 1
print(k)

parameters = c(c.m = c.m, c.s = c.s, beta.m = beta.m, beta.s = beta.s, n = n,
               p.m = p.m, p.s = p.s, epsilon = epsilon, k=as.numeric(k))
V.m <- inoculum.factor*initial.titre; T.m <- 1; T.s <- 1 # equations are proportions
state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s, 0)
state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                 "V.h","T.s","I.s")
names(state) <- state.names
                 
dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
## plot(dat[,1],1 * (dat[,"V.m"] + dat[,"V.h"]))
titre.disseminated <- rowSums(dat[,c("V.m","V.h")])

## Get model output for failed infections
titre.disseminated <- rowSums(dat[data.days,c("V.m","V.h")])
titre.data <- floor(10^fu$log10.mean.titre[model.days] + 0.5)
titre.failure <- (V.m*exp(-c.m*times))[data.days]
prop.failure <- pmax(fu$Detection.rate[model.days] - prop.disseminated,0)
denom <- prop.failure + prop.disseminated
titre.total <- (titre.disseminated*prop.disseminated + titre.failure * prop.failure) / denom

## generate plot
tiff(paste0("../figures/model_output_fixed_EIP_Fu_inoculum_x",inoculum.factor,".tif"),
            res=600, compression="lzw",height=600*6,width=600*6)
par(mar=c(5.1,4.1,4.1,2.1))
layout(mat = matrix(c(1:6,rep(7,3)), nrow = 3,byrow=TRUE))
plot(dat[,1]/24, dat[,"V.m"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'V.m', las = 1)
plot(dat[,1]/24, dat[,"T.m"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'T.m', las = 1)
plot(dat[,1]/24, dat[,"I.mn"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'I.mn', las = 1)
plot(dat[,1]/24, dat[,"V.h"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'V.h', las = 1)
plot(dat[,1]/24, dat[,"T.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'T.s', las = 1)
plot(dat[,1]/24, dat[,"I.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'I.s0', las = 1)
plot(dat[,1]/24,dat[,"V.m"] + dat[,"V.h"],lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^0.75,10^fu$log10.mean.titre, dat[,"V.m"] + dat[,"V.h"])))
points(fu$Time.pi..hours/24,10^fu$log10.mean.titre)
lines(dat[,1]/24, dat[,"V.m"],col="red",lwd=2)
lines(dat[,1]/24, dat[,"V.h"],col="green",lwd=2)
abline(h=10^0.75,lty="dashed")
dev.off()

###  Fit with mle2 - do this jointly with secondary stage???
## fixed parameters
times = seq(from = 0, to = 14*24, by = 1)
pars.init <- c(beta.m=0,p.m=0,c.m=0,epsilon=0)
beta.s <- exp(optim.out.intrathoracic$par["beta"])/24
p.s <- exp(optim.out.intrathoracic$par["p.s"])/24
c.s <- p.s / geometric.mean(tail(fu.intrathoracic$titre))

NLL <- function(pars, times) {
    ## n <- floor(as.numeric(exp(pars[5]))+0.5)
    beta.m <- as.numeric(exp(pars[1]))
    p.m <- as.numeric(exp(pars[2]))
    c.m <- as.numeric(exp(pars[3]))
    epsilon <- as.numeric(exp(pars[4]))
    k <- geometric.mean(tail(10^fu$log10.mean.titre))/geometric.mean(tail(fu.intrathoracic$titre)) -
        p.m*(1-exp(-beta.m*initial.titre))/p.s
    if (k<0) {return(Inf)}
    parameters = c(c.m = c.m, c.s = c.s, beta.m = beta.m, beta.s = beta.s, n = n,
                   p.m = p.m, p.s = p.s, epsilon = epsilon, k=as.numeric(k))
    V.m <- initial.titre; T.m <- 1; T.s <- 1 # equations are proportions
    state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s, 0)
    state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                 "V.h","T.s","I.s")
    names(state) <- state.names
    dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
    data.days <- which(dat[,"time"] %in% fu$Time.pi..hours.)
    model.days <- which(fu$Time.pi..hours. %in% dat[,"time"])
    titre.disseminated <- rowSums(dat[data.days,c("V.m","V.h")])
    titre.data <- floor(10^fu$log10.mean.titre[model.days] + 0.5)
    titre.failure <- (V.m*exp(-c.m*times))[data.days]
    prop.failure <- pmax(fu$Detection.rate[model.days] - prop.disseminated,0)
    denom <- prop.failure + prop.disseminated
    titre.total <- (titre.disseminated*prop.disseminated + titre.failure * prop.failure) / denom
    ## return(-sum(dpois(titre.data,titre.disseminated,log=TRUE)))
    return(-sum(dpois(titre.data,titre.total,log=TRUE)))
}

optim.out <- mle2(pars.init,NLL,times=times)
save(optim.out,file="optim_out.RData")


### Boneyard
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
##     data.days <- which(dat[,"time"] %in% fu$Time.pi..hours.)
##     model.days <- which(fu$Time.pi..hours. %in% dat[,"time"])
##     logtitre.model <- log((dat[,"V.m"] + dat[,"V.h"])[data.days],10)
##     logtitre.data <- fu$log10.mean.titre[model.days]
##     return(sum(abs(logtitre.model - logtitre.data)))
## }
