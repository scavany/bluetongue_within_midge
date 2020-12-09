## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspace
rm(list = ls())

## install necessary packages
if(!require(deSolve)){install.packages('deSolve'); library(deSolve)}
if(!require(BayesianTools)){install.packages('BayesianTools'); library(BayesianTools)}
if(!require(viridis)){install.packages('viridis'); library(viridis)}

## load within-vector model of BTV infection and other functions
source("./general_functions.R")
source("./btv_culicoides_functions.R")
source("./within_midge_staged_fn.R")

### Parameters - assumed
Dt <- 1/24
tmax <- 61
tvec <- seq(0,tmax,Dt)
Tmean <- 21
Tamp <- 11
k <- sqrt(1e7)
epsilon <- 2
n <- 2
Tvec <- temperature.curve(tvec,Tmean,Tamp)

###Parameters - literature derived
vector.competence.10 <- sin(18.2343*pi/180)^2
vector.competence.16 <- sin(20.6322*pi/180)^2

## set-up model to get viral abundance curves
load("./optim_out.RData",verbose=T)
## pars <- as.numeric(MAP(out)$parametersMAP) # if mcmc
pars <- as.numeric(optim.out$par)
n <- as.integer(exp(pars[1])+0.5)
beta <- exp(pars[2])
p.m <- exp(pars[3])
p.s <- p.m*exp(pars[4])
c <- (p.s + p.m) / k
parameters = c(c.m = c, c.h = c, beta.m = beta, beta.s = beta, n = n,
               p.m = p.m, p.s = p.s, epsilon = 2)
V.m <- 1; T.m <- 1; T.s <- 1 # equations are proportions
state = c(V.m, T.m, rep(0,parameters["n"]+1),0, T.s, rep(0,parameters["n"]+1))
state.names <- c("V.m","T.m",paste0("I.m",0:(parameters["n"]-1)),"I.mn",
                 "V.h","T.s",paste0("I.s",0:(parameters["n"]-1)),"I.sn")
names(state) <- state.names

## run model to get viral abundance curves
incubation.period.fn.10 <- function(t,mean=Tmean,amp=Tamp) {
    T <- temperature.curve(time=t,mean=mean,amp=amp)
    return(incubation.period.10(T))
}
dat.10 = ode(y = state, times = tvec, func = wv.BTV.varT, parms = parameters,
          incubation.period.fn = incubation.period.fn.10)
incubation.period.fn.16 <- function(t,mean=Tmean,amp=Tamp) {
    T <- temperature.curve(time=t,mean=mean,amp=amp)
    return(incubation.period.16(T))
}
dat.16 = ode(y = state, times = tvec, func = wv.BTV.varT, parms = parameters,
          incubation.period.fn = incubation.period.fn.16)

### plot the ODE
## Mellor "data"
mellor <- read.csv("../data/Mellor_points.csv",header=F)
names(mellor) <- c("days","titre")
mellor$days <- floor(mellor$days+0.5)
initial.titre <- 10^mellor$titre[1]

## plot 10
t.max.plot.10 <- which.min(abs(dat.10[,1] - 13.1))
tiff("../figures/model_output_temperature_10.tif",res=600,
     compression="lzw",height=600*6,width=600*6)
par(mar=c(5.1,4.1,4.1,2.1))
layout(mat = matrix(c(1:8,rep(9,4)), nrow = 3,byrow=TRUE))
plot(dat.10[1:t.max.plot.10,1], dat.10[1:t.max.plot.10,"V.m"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'V.m', las = 1)
plot(dat.10[1:t.max.plot.10,1], dat.10[1:t.max.plot.10,"T.m"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'T.m', las = 1)
plot(dat.10[1:t.max.plot.10,1], dat.10[1:t.max.plot.10,"I.m0"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'I.m0', las = 1)
plot(dat.10[1:t.max.plot.10,1], dat.10[1:t.max.plot.10,"I.mn"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'I.mn', las = 1)
plot(dat.10[1:t.max.plot.10,1], dat.10[1:t.max.plot.10,"V.h"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'V.h', las = 1)
plot(dat.10[1:t.max.plot.10,1], dat.10[1:t.max.plot.10,"T.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'T.s', las = 1)
plot(dat.10[1:t.max.plot.10,1], dat.10[1:t.max.plot.10,"I.s0"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'I.s0', las = 1)
plot(dat.10[1:t.max.plot.10,1], dat.10[1:t.max.plot.10,"I.sn"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'I.sn', las = 1)
plot(dat.10[1:t.max.plot.10,1],
     initial.titre*(dat.10[1:t.max.plot.10,"V.m"] + dat.10[1:t.max.plot.10,"V.h"]),
     lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^mellor$titre,
                    initial.titre*(dat.10[1:t.max.plot.10,"V.m"] + dat.10[1:t.max.plot.10,"V.h"]))))
points(mellor$day,10^mellor$titre)
dev.off()

## plot 16
t.max.plot.16 <- which.min(abs(dat.16[,1] - 13.1))
tiff("../figures/model_output_temperature_16.tif",res=600,
     compression="lzw",height=600*6,width=600*6)
par(mar=c(5.1,4.1,4.1,2.1))
layout(mat = matrix(c(1:8,rep(9,4)), nrow = 3,byrow=TRUE))
plot(dat.16[1:t.max.plot.16,1], dat.16[1:t.max.plot.16,"V.m"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'V.m', las = 1)
plot(dat.16[1:t.max.plot.16,1], dat.16[1:t.max.plot.16,"T.m"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'T.m', las = 1)
plot(dat.16[1:t.max.plot.16,1], dat.16[1:t.max.plot.16,"I.m0"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'I.m0', las = 1)
plot(dat.16[1:t.max.plot.16,1], dat.16[1:t.max.plot.16,"I.mn"], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'I.mn', las = 1)
plot(dat.16[1:t.max.plot.16,1], dat.16[1:t.max.plot.16,"V.h"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'V.h', las = 1)
plot(dat.16[1:t.max.plot.16,1], dat.16[1:t.max.plot.16,"T.s"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'T.s', las = 1)
plot(dat.16[1:t.max.plot.16,1], dat.16[1:t.max.plot.16,"I.s0"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'I.s0', las = 1)
plot(dat.16[1:t.max.plot.16,1], dat.16[1:t.max.plot.16,"I.sn"], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'I.sn', las = 1)
plot(dat.16[1:t.max.plot.16,1],
     initial.titre*(dat.16[1:t.max.plot.16,"V.m"] + dat.16[1:t.max.plot.16,"V.h"]),
     lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^mellor$titre,
                    initial.titre*(dat.16[1:t.max.plot.16,"V.m"] + dat.16[1:t.max.plot.16,"V.h"]))))
points(mellor$day,10^mellor$titre)
dev.off()

### Plot the original figure in grant but with new reassortment method
survival.probs <- vector(mode="numeric",length=length(Tvec))
incubation.probs <- survival.probs; reassortment.probs <- incubation.probs
reassortment.proportionality <- 1e-8
for (i in 1:length(Tvec)) {
    print(i)
    survival.probs[i] <- 1-rate.prob(hazard.lowRH,Tvec[1:i],Dt)
    incubation.probs[i] <- period.prob(incubation.period.10,Tvec[1:i],Dt)
    ## placeholder for reassortment - assume rate proportional to product of abundances
    reassortment.probs[i] <- 1 - exp(- ifelse(i>1,
                                              integrate(reassortment.proportionality
                                                        *rowSums(dat.10[1:i,c("V.m","V.h")])
                                                        *rowSums(dat.16[1:i,c("V.m","V.h")]),
                                                        Dt),
                                              integrate(reassortment.proportionality
                                                        *sum(dat.10[1:i,c("V.m","V.h")])
                                                        *sum(dat.16[1:i,c("V.m","V.h")]),
                                                        Dt)))
}

##plot
tiff("../figures/rates_and_probabilities_new.tif",res=600,
     compression="lzw",height=600*6,width=600*6)
par(mfrow=c(5,1),mar=c(2.1,4.1,1.1,4.1),oma=c(3,0,3,0))
plot(tvec, survival.probs, type='l', ylim=c(0,1),xlab="",ylab="Survival prob.")
par(new=TRUE)
plot(tvec,hazard.lowRH(Tvec),type='l',col="gray",yaxt="n",ylab="",xaxt="n",xlab="")
axis(4,labels=T,col="gray",col.ticks="gray",col.axis="gray")
mtext("Mortality rate",4,3,col="gray",cex=2/3)
plot(tvec, incubation.probs, type='l', ylim=c(0,1),xlab="",ylab="Incubation prob.")
par(new=TRUE)
plot(tvec,incubation.period.10(Tvec),type='l',col="gray",yaxt="n",ylab="",xaxt="n",xlab="")
axis(4,labels=T,col="gray",col.ticks="gray",col.axis="gray")
mtext("Incubation period",4,3,col="gray",cex=2/3)
plot(tvec, reassortment.probs, type='l', ylim=c(0,1),xlab="",ylab="Reassortment prob.")
par(new=TRUE)
reassortment.rates <- reassortment.proportionality*rowSums(dat.10[,c("V.m","V.h")])*rowSums(dat.16[,c("V.m","V.h")])
plot(tvec,reassortment.rates,type='l',col="gray",yaxt="n",ylab="",xaxt="n",xlab="")
axis(4,labels=T,col="gray",col.ticks="gray",col.axis="gray")
mtext("Reassortment rate",4,3,col="gray",cex=2/3)
plot(tvec, biting.rate(Tvec), type='l',xlab="",ylab="Expected bites / day")
reassortant.bites <- biting.rate(Tvec)*survival.probs*incubation.probs*reassortment.probs*vector.competence.10
plot(tvec, reassortant.bites,
     type='l',xlab="",ylab="Infectious bites")
mtext("with reassortant / day",2,2,cex=2/3)
polygon(c(tvec,rev(tvec)),c(reassortant.bites,rep(0,length(reassortant.bites))),col="red",border=F)
dev.off()
