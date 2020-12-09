library(data.table)
## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts')

## clear existing workspace
rm(list = ls())

## install necessary packages
if(!require(deSolve)){install.packages('deSolve'); library(deSolve)}

## specify within-vector model of BTV infection
wv.BTV = function(t, state, parameters)
{
    with(as.list(c(state, parameters)),{
        dV = p*I2 - c*V
        dT = -beta*T*V
        dI1 = beta*T*V - I1/eclipse
        dI2 = I1/eclipse
        list(c(dV, dT, dI1, dI2))
    })
}

## run model
V0 <- 100
Vinf <- V0*10^3.5
eclipse <- 1
times = seq(from = 0, to = 13, by = 0.01)
## parms = c(beta = 1e-5, c = 1, p = 1, sigma = 0.1) 
parms.list <- expand.grid(list(beta=10^seq(-5,0),c=10^seq(-5,0),p=10^seq(-5,0),sigma=10^seq(-5,0)))
## state = c(V=V0,T=parms["c"]*Vinf/parms["p"],I=0)
state.list = data.frame(V = rep(V0,nrow(parms.list)), T = Vinf*parms.list$c/parms.list$p,
                        I = rep(0,nrow(parms.list)))

time.vec <- rep(NA,nrow(parms.list))
min.vec <-  rep(NA,nrow(parms.list))
min.time.vec <- rep(NA,nrow(parms.list))
for (i in seq(1,nrow(parms.list))) {
    print(i)
    parms <- unlist(parms.list[i,])
    state <- unlist(state.list[i,])
    dat = ode(y = state, times = times, func = wv.BTV, parms = parms)
    min.vec[i] <- min(dat[,"V"])
    min.time.vec[i] <- dat[which.min(dat[,"V"]),"time"]
    time.vec[i] <- dat[min(which(dat[,"V"] > 1e5)),"time"]
}

parm.options <- parms.list[min.vec<87,]
save(parm.options,file="parm_options.RData")

## generate plot
load("parm_options.RData")
row <- sample(nrow(parm.options),1)
parms <- unlist(parm.options[row,])
T0 <- as.numeric(Vinf*parms["c"]/parms["p"])
state <- c(V=V0,T=T0 ,I1=0,I2=0)
dat = ode(y = state, times = times, func = wv.BTV, parms = parms)
par(mfrow=c(2,2))
plot(log(V)/log(10)~time,data=dat, type = 'l', lwd = 2, bty = 'n', yaxs="i",xaxs="i",
     xlab = 'Time (Days)', ylab = 'Viral load', las = 1)
plot(T~time,data=dat, type = 'l', lwd = 2, bty = 'n', yaxs="i",xaxs="i",
     xlab = 'Time (Days)', ylab = 'Target cells', las = 1)
plot(I1~time,data=dat, type = 'l', lwd = 2, bty = 'n', yaxs="i",xaxs="i",
     xlab = 'Time (Days)', ylab = 'Eclipse phase cells', las = 1)
plot(I2~time,data=dat, type = 'l', lwd = 2, bty = 'n', yaxs="i",xaxs="i",
     xlab = 'Time (Days)', ylab = 'Infected cells', las = 1)



## coinfection model
times1 <- times
times2 <- times + 5
virus1 <- ode(y = state, times = times1, func = wv.BTV, parms = parms)
virus2 <- ode(y = state, times = times2, func = wv.BTV, parms = parms)

