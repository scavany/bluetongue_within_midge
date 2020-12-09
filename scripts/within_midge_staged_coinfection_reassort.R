## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspace
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(deSolve, BayesianTools, data.table, viridis,mc2d)

## load within-vector model of BTV infection
source("./general_functions.R")
source("./btv_culicoides_functions.R")
source("./within_midge_staged_fn.R")
load("./optim_out.RData")

## Mellor "data"
mellor <- read.csv("../data/Mellor_points.csv",header=F)
names(mellor) <- c("days","titre")
mellor$days <- floor(mellor$days+0.5)
initial.titre <- 10^mellor$titre[1]

## Samal data
samal <- fread("../data/samal_data.csv")

## el hussein data
el.hussein <- fread("../data/elHussein_data.csv")
##plot elHussein data
tiff("../figures/hussein_data.tif",res=600,
     compression="lzw",height=400*9,width=400*16)
par(mfrow=c(2,1),mar=c(3.1,4.1,1.1,2.1),
    oma=c(3,0,3,0))
cols=viridis(10)
plot(-1,-1,xlim=c(0,15),ylim=c(0,1),xlab="",
     ylab="Proportion reassortant",bty="n")
points(el.hussein[gap==0]$time, el.hussein[gap==0]$prop_reassortant,col=cols[1],pch=19)
points(el.hussein[gap==1]$time, el.hussein[gap==1]$prop_reassortant,col=cols[2],pch=19)
points(el.hussein[gap==3]$time, el.hussein[gap==3]$prop_reassortant,col=cols[4],pch=19)
points(el.hussein[gap==5]$time, el.hussein[gap==5]$prop_reassortant,col=cols[6],pch=19)
points(el.hussein[gap==7]$time, el.hussein[gap==7]$prop_reassortant,col=cols[8],pch=19)
points(el.hussein[gap==9]$time, el.hussein[gap==9]$prop_reassortant,col=cols[10],pch=19)
legend("topleft",legend=paste(c(0,seq(1,9,2)),"day gap"),bty="n",col=cols[1+c(0,seq(1,9,2))],
       pch=19)
cols=viridis(7)
plot(-1,-1,xlim=c(0,15),ylim=c(0,1),xlab="",
     ylab="Proportion reassortant",bty="n")
points(el.hussein[time==5]$gap, el.hussein[time==5]$prop_reassortant,col=cols[1],pch=19)
points(el.hussein[time==7]$gap, el.hussein[time==7]$prop_reassortant,col=cols[3],pch=19)
points(el.hussein[time==9]$gap, el.hussein[time==9]$prop_reassortant,col=cols[5],pch=19)
points(el.hussein[time==10]$gap, el.hussein[time==10]$prop_reassortant,col=cols[6],pch=19)
points(el.hussein[time==11]$gap, el.hussein[time==11]$prop_reassortant,col=cols[7],pch=19)
legend("topright",legend=paste(c(5,7,9,10,11),"days since first infection"),bty="n",
       col=cols[-4+c(5,7,9,10,11)],
       pch=19)
mtext("Time (Days)",side=1,line=3,cex=3/4)
dev.off()

## Load mcmc output
load("mcmc_reassort_out.RData",verbose=TRUE)
samples = getSample(out,thin=10,start=1000)
growth.ratio <- colMeans(samples)[1]
withlike <- colMeans(samples)[2]

## fixed parameters
pars <- as.numeric(optim.out$par)
n <- as.integer(exp(pars[1])+0.5)
beta <- exp(pars[2])
p.m <- exp(pars[3])
p.s <- p.m*exp(pars[4])
k <- sqrt(1e7)
epsilon <- 3
eip <- 10
tau <- eip - epsilon
## run model
times = seq(from = 0, to = 20, by = 0.1)
c <- (p.s + p.m) / k
parameters = c(c.m = c, c.s = c, beta.m = beta, beta.s = beta, n = n,
               p.m = p.m, p.s = p.s, epsilon = epsilon, tau = tau)
V.ma <- 1; V.mb <- 1; T.m <- 1; T.s <- 1 # equations are proportions
state = c(V.ma, V.mb,
          T.m, rep(0,parameters["n"]^2 + 4*parameters["n"] + 1),
          rep(0,3),
          T.s, rep(0,3*parameters["n"]*(parameters["n"] + 3)))
state.names <- c("V.ma","V.mb","T.m",
                 paste0("I.ma",0:(parameters["n"]-1)),"I.man",
                 paste0("I.mb",0:(parameters["n"]-1)),"I.mbn", "I.mab00",
                 paste0("I.mab",1:(parameters["n"]-1),0),
                 paste0("I.mab",0,1:(parameters["n"]-1)),
                 paste0("I.mab",rep(1:(parameters["n"]-1),parameters["n"]-1),
                        rep(1:(parameters["n"]-1),each=parameters["n"]-1)),
                 paste0("I.mabn",1:(parameters["n"]-1)),
                 paste0("I.mab",1:(parameters["n"]-1),"n"),
                 "I.mabnn",
                 "V.sa","V.sb","V.sr","T.s",
                 paste0("I.sa",0:(parameters["n"]-1)),"I.san",
                 paste0("I.sb",0:(parameters["n"]-1)),"I.sbn",
                 paste0("I.sr",0:(parameters["n"]-1)),"I.srn",
                 "I.sab00",
                 paste0("I.sab",1:(parameters["n"]-1),0),
                 paste0("I.sab",0,1:(parameters["n"]-1)),
                 paste0("I.sab",rep(1:(parameters["n"]-1),parameters["n"]-1),
                        rep(1:(parameters["n"]-1),each=parameters["n"]-1)),
                 paste0("I.sabn",1:(parameters["n"]-1)),
                 paste0("I.sab",1:(parameters["n"]-1),"n"),
                 "I.sabnn",
                 "I.sar00",
                 paste0("I.sar",1:(parameters["n"]-1),0),
                 paste0("I.sar",0,1:(parameters["n"]-1)),
                 paste0("I.sar",rep(1:(parameters["n"]-1),parameters["n"]-1),
                        rep(1:(parameters["n"]-1),each=parameters["n"]-1)),
                 paste0("I.sarn",1:(parameters["n"]-1)),
                 paste0("I.sar",1:(parameters["n"]-1),"n"),
                 "I.sarnn",
                 "I.sbr00",
                 paste0("I.sbr",1:(parameters["n"]-1),0),
                 paste0("I.sbr",0,1:(parameters["n"]-1)),
                 paste0("I.sbr",rep(1:(parameters["n"]-1),parameters["n"]-1),
                        rep(1:(parameters["n"]-1),each=parameters["n"]-1)),
                 paste0("I.sbrn",1:(parameters["n"]-1)),
                 paste0("I.sbr",1:(parameters["n"]-1),"n"),
                 "I.sbrnn")
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
                     func = wv.BTV.coinfection.reassort, parms = parameters,
                     prod.fn.a=prod.wrapper.a,
                     prod.fn.b=prod.wrapper.b,
                     prod.fn.r=prod.wrapper.r))

## generate plot
tiff("../figures/model_output_fixed_EIP_coinfection_reassort.tif",res=600,
     compression="lzw",height=600*7,width=600*7)
par(mar=c(2.1,4.1,1.1,2.1),
    oma=c(3,0,3,0))
layout(mat = matrix(c(1,1,2,2), nrow = 2,byrow=TRUE))
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
     xlab = '', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^mellor$titre, 0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]))))
lines(dat[,1],0.5*initial.titre*(dat[,"V.mb"] + dat[,"V.sb"]),col="gray",lty="dashed",lwd=2)
lines(dat[,1],0.5*initial.titre*(dat[,"V.sr"]),col="red",lwd=2)
lines(dat[,1],0.5*initial.titre*(rowSums(dat[,grepl("V",names(dat))])),col="green",lwd=1)
points(mellor$day,10^mellor$titre)
legend("topleft",bty="n",legend=c("Strain A","Strain B","Reassortants", "Total"),
       lty=c("solid","dashed","solid","solid"),col=c("black","gray","red","green"),
       lwd=c(2,2,2,1))
plot(dat[,1],dat[,"V.sr"]/rowSums(dat[,grepl("V",names(dat))]),lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Proportion reassortant', las = 1,type='l',ylim=c(0,1))
mtext("Time (Days)",side=1,line=3,cex=3/4)
points(samal$time,samal$prop_reassortant,pch=1)
points(el.hussein[gap==0]$time, el.hussein[gap==0]$prop_reassortant,pch=2)
legend("topleft",bty="n",legend=c("Samal, 1987","el Hussein, 1989"),
       pch=c(1,2))
dev.off()


### With temperature dependence
## extra parameters
Dt <- 1/24
tmax <- 20
tvec <- seq(0,tmax,Dt)
Tmean <- 21
Tamp <- 11
incubation.period.fn.10 <- function(t,mean=Tmean,amp=Tamp) {
    T <- temperature.curve(time=t,mean=mean,amp=amp)
    return(incubation.period.10(T))
}

## reassortment parameters
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
                     func = wv.BTV.coinfection.reassort.varT, parms = parameters,
                     incubation.period.fn = incubation.period.fn.10,
                     prod.fn.a=prod.wrapper.a,
                     prod.fn.b=prod.wrapper.b,
                     prod.fn.r=prod.wrapper.r))

tiff("../figures/model_output_coinfection_reassort_varT.tif",res=600,
     compression="lzw",height=600*7,width=600*7)
par(mar=c(2.1,4.1,1.1,2.1),
    oma=c(3,0,3,0))
layout(mat = matrix(c(1,1,2,2), nrow = 2,byrow=TRUE))
plot(dat[,1],0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]),lwd = 2, bty = 'n', 
     xlab = '', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^mellor$titre, 0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]))))
lines(dat[,1],0.5*initial.titre*(dat[,"V.mb"] + dat[,"V.sb"]),col="gray",lty="dashed",lwd=2)
lines(dat[,1],0.5*initial.titre*(dat[,"V.sr"]),col="red",lwd=2)
lines(dat[,1],0.5*initial.titre*(rowSums(dat[,grepl("V",names(dat))])),col="green",lwd=1)
points(mellor$day,10^mellor$titre)
legend("topleft",bty="n",legend=c("Strain A","Strain B","Reassortants", "Total"),
       lty=c("solid","dashed","solid","solid"),col=c("black","gray","red","green"),
       lwd=c(2,2,2,1))
plot(dat[,1],dat[,"V.sr"]/rowSums(dat[,grepl("V",names(dat))]),lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Proportion reassortant', las = 1,type='l',ylim=c(0,1))
mtext("Time (Days)",side=1,line=3,cex=3/4)
points(samal$time,samal$prop_reassortant,pch=1)
points(el.hussein[gap==0]$time, el.hussein[gap==0]$prop_reassortant,pch=2)
points(0.95*max(dat[,1]),
       mean(c(samal$prop_reassortant,el.hussein[gap==0]$prop_reassortant)),
       col="red",pch=23)
legend("topleft",bty="n",legend=c("Samal, 1987","el Hussein, 1989", "Mean"),
       pch=c(1,2,23),col=c("black","black","red"))
dev.off()


### Different times
for (second.intro in (1:5)*2-1){
    times.1 <- seq(0,second.intro,0.1)
    times.2 <- seq(second.intro,20,0.1)
    state.1 <- state
    state.1["V.mb"] <- 0
    dat.1 = data.frame(ode(y = state.1, times = times.1,
                           func = wv.BTV.coinfection.reassort.varT, parms = parameters,
                           incubation.period.fn = incubation.period.fn.10,
                           prod.fn.a=prod.wrapper.a,
                           prod.fn.b=prod.wrapper.b,
                           prod.fn.r=prod.wrapper.r))
    state.2 <- as.numeric(dat.1[nrow(dat.1),-1])
    names(state.2) <- names(state)
    state.2["V.mb"] <- 1
    dat.2 = data.frame(ode(y = state.2, times = times.2,
                           func = wv.BTV.coinfection.reassort.varT, parms = parameters,
                           incubation.period.fn = incubation.period.fn.10,
                           prod.fn.a=prod.wrapper.a,
                           prod.fn.b=prod.wrapper.b,
                           prod.fn.r=prod.wrapper.r))
    dat <- rbind(dat.1,dat.2)
    tiff(paste0("../figures/model_output_coinfection_reassort_varT_sep",second.intro,
                "_",withlike,"wl.tif"),
         res=600,compression="lzw",height=600*7,width=600*7)
    par(mar=c(2.1,4.1,1.1,2.1),
        oma=c(3,0,3,0))
    layout(mat = matrix(c(1,1,2,2), nrow = 2,byrow=TRUE))
    plot(dat[,1],0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]),lwd = 2, bty = 'n', 
         xlab = '', ylab = 'Viral load', las = 1,log="y",type='l',
         ylim = range(c(10^mellor$titre, 0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]))))
    lines(dat[,1],0.5*initial.titre*(dat[,"V.mb"] + dat[,"V.sb"]),col="gray",lty="dashed",lwd=2)
    lines(dat[,1],0.5*initial.titre*(dat[,"V.sr"]),col="red",lwd=2)
    lines(dat[,1],0.5*initial.titre*(rowSums(dat[,grepl("V",names(dat))])),col="green",lwd=1)
    points(mellor$day,10^mellor$titre)
    legend("topleft",bty="n",legend=c("Strain A","Strain B","Reassortants", "Total"),
           lty=c("solid","dashed","solid","solid"),col=c("black","gray","red","green"),
           lwd=c(2,2,2,1))
    plot(dat[,1],dat[,"V.sr"]/rowSums(dat[,grepl("V",names(dat))]),lwd = 2, bty = 'n', 
         xlab = 'Time (Days)', ylab = 'Proportion reassortant', las = 1,type='l',ylim=c(0,1))
    mtext("Time (Days)",side=1,line=3,cex=3/4)
    points(el.hussein[gap==second.intro]$time, el.hussein[gap==second.intro]$prop_reassortant,pch=2)
    legend("topleft",bty="n",legend=c("el Hussein, 1989"),
           pch=c(2),col=c("black"))
    dev.off()
}
