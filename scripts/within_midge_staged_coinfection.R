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
epsilon <- 3
eip <- 10
tau <- eip - epsilon
n <- 3
## run model
times = seq(from = 0, to = 13, by = 0.1)
beta <- 1.42
p.m <- 2.43
p.s <- p.m*2950 ## p.s > p.m
c <- (p.s + p.m) / k
growth.ratio <- 100^3
withlike <- 0.75
parameters = c(c.m = c, c.s = c, beta.m = beta, beta.s = beta, n = n,
               p.m = p.m, p.s = p.s, epsilon = epsilon, tau = tau)
V.ma <- 1; V.mb <- 1; T.m <- 1; T.s <- 1 # equations are proportions
state = c(V.ma, V.mb,
          T.m, rep(0,parameters["n"]^2 + 4*parameters["n"] + 1),
          rep(0,3),
          T.s, rep(0,parameters["n"]^2 + 4*parameters["n"] + 1))
state.names <- c("V.ma","V.mb","T.m",
                 paste0("I.ma",0:(parameters["n"]-1)),"I.man",
                 paste0("I.mb",0:(parameters["n"]-1)),"I.mbn", "I.mab00",
                 paste0("I.mab",1:(parameters["n"]-1),0),
                 paste0("I.mab",0,1:(parameters["n"]-1)),
                 paste0("I.mab",rep(1:(parameters["n"]-1),parameters["n"]-1),
                        rep(1:(parameters["n"]-1),each=parameters["n"]-1)),
                 paste0("I.mabn",1:(parameters["n"]-1)),
                 paste0("I.mab",1:(parameters["n"]-1),"n"),
                 "I.mcn",
                 "V.sa","V.sb","V.sr","T.s",
                 paste0("I.sa",0:(parameters["n"]-1)),"I.san",
                 paste0("I.sb",0:(parameters["n"]-1)),"I.sbn", "I.sab00",
                 paste0("I.sab",1:(parameters["n"]-1),0),
                 paste0("I.sab",0,1:(parameters["n"]-1)),
                 paste0("I.sab",rep(1:(parameters["n"]-1),parameters["n"]-1),
                        rep(1:(parameters["n"]-1),each=parameters["n"]-1)),
                 paste0("I.sabn",1:(parameters["n"]-1)),
                 paste0("I.sab",1:(parameters["n"]-1),"n"),
                 "I.scn")
names(state) <- state.names
prod.wrapper.a <- function(i,j,n,r = growth.ratio, wl = withlike){
    production.fn.a(i,j,n,r,wl)
}
prod.wrapper.b <- function(i,j,n,r = growth.ratio, wl = withlike){
    production.fn.b(i,j,n,r,wl)
}
prod.wrapper.r <- function(i,j,n,r = growth.ratio, wl = withlike){
    production.fn.r(i,j,n,r,wl)
}           
dat = data.frame(ode(y = state, times = times,
                     func = wv.BTV.coinfection, parms = parameters,
                     prod.fn.a=prod.wrapper.a,
                     prod.fn.b=prod.wrapper.b,
                     prod.fn.r=prod.wrapper.r))

## generate plot
tiff("../figures/model_output_fixed_EIP_coinfection_new.tif",res=600,
     compression="lzw",height=600*6,width=600*12)
## par(mar=c(5.1,4.1,4.1,2.1))
## layout(mat = matrix(c(1:8,rep(9,4)), nrow = 3,byrow=TRUE))
## plot(dat[,1], dat[,"V.m"], type = 'l', lwd = 2, bty = 'n', 
##      xlab = 'Time (Days)', ylab = 'V.m', las = 1)
## plot(dat[,1], dat[,"T.m"], type = 'l', lwd = 2, bty = 'n',
##      xlab = 'Time (Days)', ylab = 'T.m', las = 1)
## plot(dat[,1], dat[,"I.m0"], type = 'l', lwd = 2, bty = 'n',
##      xlab = 'Time (Days)', ylab = 'I.m0', las = 1)
## plot(dat[,1], dat[,"I.mn"], type = 'l', lwd = 2, bty = 'n',
##      xlab = 'Time (Days)', ylab = 'I.mn', las = 1)
## plot(dat[,1], dat[,"V.h"], type = 'l', lwd = 2, bty = 'n', 
##      xlab = 'Time (Days)', ylab = 'V.h', las = 1)
## plot(dat[,1], dat[,"T.s"], type = 'l', lwd = 2, bty = 'n', 
##      xlab = 'Time (Days)', ylab = 'T.s', las = 1)
## plot(dat[,1], dat[,"I.s0"], type = 'l', lwd = 2, bty = 'n', 
##      xlab = 'Time (Days)', ylab = 'I.s0', las = 1)
## plot(dat[,1], dat[,"I.sn"], type = 'l', lwd = 2, bty = 'n', 
##      xlab = 'Time (Days)', ylab = 'I.sn', las = 1)
plot(dat[,1],0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]),lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^mellor$titre, 0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]))))
lines(dat[,1],0.5*initial.titre*(dat[,"V.mb"] + dat[,"V.sb"]),col="gray",lty="dashed",lwd=2)
lines(dat[,1],0.5*initial.titre*(dat[,"V.sr"]),col="red",lwd=2)
lines(dat[,1],0.5*initial.titre*(rowSums(dat[,grepl("V",names(dat))])),col="green",lwd=1)
points(mellor$day,10^mellor$titre)
legend("topleft",bty="n",legend=c("Strain A","Strain B","Reassortants", "Total"),
       lty=c("solid","dashed","solid","solid"),col=c("black","gray","red","green"),
       lwd=c(2,2,2,1))
dev.off()





###  Fit with optim - same as before at the moment - needs work!!!!!
## fixed parameters
times = seq(from = 0, to = 13, by = 0.1)
pars.init <- c(n=log(3),beta=log(0.92),p.m=log(5.7),p.s.mult=log(890))

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
    parameters = c(c.m = c, c.s = c, beta.m = beta, beta.s = beta, n = n,
                   p.m = p.m, p.s = p.s, epsilon = epsilon, tau = tau)
    V.ma <- 1; V.mb <- 1; T.m <- 1; T.s <- 1 # equations are proportions
    state = c(V.ma, V.mb,
              T.m, rep(0,parameters["n"]^2 + 4*parameters["n"] + 1),
              rep(0,3),
              T.s, rep(0,parameters["n"]^2 + 4*parameters["n"] + 1))
    state.names <- c("V.ma","V.mb","T.m",
                     paste0("I.ma",0:(parameters["n"]-1)),"I.man",
                     paste0("I.mb",0:(parameters["n"]-1)),"I.mbn", "I.mab00",
                     paste0("I.mab",1:(parameters["n"]-1),0),
                     paste0("I.mab",0,1:(parameters["n"]-1)),
                     paste0("I.mab",rep(1:(parameters["n"]-1),parameters["n"]-1),
                            rep(1:(parameters["n"]-1),each=parameters["n"]-1)),
                     paste0("I.mabn",1:(parameters["n"]-1)),
                     paste0("I.mab",1:(parameters["n"]-1),"n"),
                     "I.mcn",
                     "V.sa","V.sb","V.sr","T.s",
                     paste0("I.sa",0:(parameters["n"]-1)),"I.san",
                     paste0("I.sb",0:(parameters["n"]-1)),"I.sbn", "I.sab00",
                     paste0("I.sab",1:(parameters["n"]-1),0),
                     paste0("I.sab",0,1:(parameters["n"]-1)),
                     paste0("I.sab",rep(1:(parameters["n"]-1),parameters["n"]-1),
                            rep(1:(parameters["n"]-1),each=parameters["n"]-1)),
                     paste0("I.sabn",1:(parameters["n"]-1)),
                     paste0("I.sab",1:(parameters["n"]-1),"n"),
                     "I.scn")
    names(state) <- state.names   
    dat = as.data.frame(ode(y = state, times = times, func = wv.BTV.coinfection, parms = parameters))
    data.days <- which(dat[,"time"] %in% mellor$days)
    logtitre.model <- log(initial.titre*(rowSums(dat[,grepl("V",names(dat))]))[data.days],10)
    logtitre.data <- mellor$titre
    return(sum(abs(logtitre.model - logtitre.data)))
}

optim.out <- optim(pars.init,optim.fn,times=times)
save(optim.out,file="optim_out_coinfection.RData")
