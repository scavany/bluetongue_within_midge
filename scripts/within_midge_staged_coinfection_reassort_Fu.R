## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspace
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(deSolve, BayesianTools, data.table, viridis,mc2d,
       grDevices,psych,extraDistr)

## load within-vector model of BTV infection
source("./general_functions.R")
source("./btv_culicoides_functions.R")
source("./within_midge_staged_fn.R")
load("./optim_out.RData",verbose=TRUE)
## load("./optim_out_intrathoracic.RData",verbose=TRUE)

## Mellor "data"
## mellor <- read.csv("../data/Mellor_points.csv",header=F)
## names(mellor) <- c("days","titre")
## mellor$days <- floor(mellor$days+0.5)
## initial.titre <- 10^mellor$titre[1]

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
fu$log10.mean.titre <- log(log(2) * 10^fu$log10.mean.titre,10) # convert to pfu
fu.hours <- fu$Time.pi..hours.
initial.titre <- 10^fu$log10.mean.titre[1]
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic$titre <- fu.intrathoracic$titre * log(2) #convert to pfu
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)

## Samal data
samal <- fread("../data/samal_data.csv")
samal$time <- samal$time*24

## el hussein data
el.hussein <- fread("../data/elHussein_data.csv")
el.hussein$time <- el.hussein$time*24
el.hussein$gap <- el.hussein$gap*24

##plot elHussein data
tiff("../figures/hussein_data.tif",res=600,
     compression="lzw",height=400*9,width=400*16)
par(mfrow=c(2,1),mar=c(3.1,4.1,1.1,2.1),
    oma=c(3,0,3,0))
cols=viridis(10)
plot(-1,-1,xlim=c(0,15*24),ylim=c(0,1),xlab="",
     ylab="Proportion reassortant",bty="n")
points(el.hussein[gap==0]$time, el.hussein[gap==0]$prop_reassortant,col=cols[1],pch=19)
points(el.hussein[gap==1*24]$time, el.hussein[gap==1*24]$prop_reassortant,col=cols[2],pch=19)
points(el.hussein[gap==3*24]$time, el.hussein[gap==3*24]$prop_reassortant,col=cols[4],pch=19)
points(el.hussein[gap==5*24]$time, el.hussein[gap==5*24]$prop_reassortant,col=cols[6],pch=19)
points(el.hussein[gap==7*24]$time, el.hussein[gap==7*24]$prop_reassortant,col=cols[8],pch=19)
points(el.hussein[gap==9*24]$time, el.hussein[gap==9*24]$prop_reassortant,col=cols[10],pch=19)
legend("topleft",legend=paste(c(0,seq(1,9,2)),"day gap"),bty="n",col=cols[1+c(0,seq(1,9,2))],
       pch=19)
cols=viridis(7)
plot(-1,-1,xlim=c(0,15*24),ylim=c(0,1),xlab="",
     ylab="Proportion reassortant",bty="n")
points(el.hussein[time==5*24]$gap, el.hussein[time==5*24]$prop_reassortant,col=cols[1],pch=19)
points(el.hussein[time==7*24]$gap, el.hussein[time==7*24]$prop_reassortant,col=cols[3],pch=19)
points(el.hussein[time==9*24]$gap, el.hussein[time==9*24]$prop_reassortant,col=cols[5],pch=19)
points(el.hussein[time==10*24]$gap, el.hussein[time==10*24]$prop_reassortant,col=cols[6],pch=19)
points(el.hussein[time==11*24]$gap, el.hussein[time==11*24]$prop_reassortant,col=cols[7],pch=19)
legend("topright",legend=paste(c(5,7,9,10,11),"days since first infection"),bty="n",
       col=cols[-4+c(5,7,9,10,11)],
       pch=19)
mtext("Time (Days)",side=1,line=3,cex=3/4)
dev.off()

## Load mcmc output
## load("mcmc_reassort_out.RData",verbose=TRUE)
## samples = getSample(out,thin=10,start=1000)
## growth.ratio <- colMeans(samples)[1] 
## withlike <- colMeans(samples)[2]
## tau <- colMeans(samples)[3]

## fixed parameters
pars <- as.numeric(optim.out$par)
n <- 1
beta.s <- as.numeric(exp(optim.out$par["lbeta.s"]))
p.s <- as.numeric(exp(optim.out$par["lp.s"]))
c.s <- as.numeric(p.s / geometric.mean(tail(fu.intrathoracic$titre)))
epsilon <- as.numeric(exp(optim.out$par["lepsilon"]))
beta.m <- as.numeric(exp(optim.out$par["lbeta.m"]))
p.m <- as.numeric(exp(optim.out$par["lp.m"]))
c.m <- as.numeric(exp(optim.out$par["lc.m"]))
k <- geometric.mean(tail(10^fu$log10.mean.titre))/geometric.mean(tail(fu.intrathoracic$titre)) -
    p.m*(1-exp(-beta.m*initial.titre))/p.s

## run model
times = seq(from = 0, to = 20*24, by = 1)
parameters = c(c.m = c.m, c.s = c.s, beta.m = beta.m, beta.s = beta.s, n = n,
               p.m = p.m, p.s = p.s, epsilon = epsilon, k = k)
V.ma <- initial.titre; V.mb <- initial.titre; T.m <- 1; T.s <- 1 # equations are proportions
state = c(V.ma, V.mb,
          T.m, rep(0,parameters["n"]^2 + 4*parameters["n"] + 1),
          rep(0,3),
          T.s, rep(0,12))
if (n>1) {
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
                     "I.sa0", "I.san",
                     "I.sb0", "I.sbn",
                     "I.sr0", "I.srn",
                     "I.sab00",
                     "I.sar00",
                     "I.sbr00",
                     "I.sabnn",
                     "I.sarnn",
                     "I.sbrnn")
} else {
    state.names <- c("V.ma","V.mb","T.m",
                     "I.ma0","I.man",
                     "I.mb0","I.mbn", "I.mab00",
                     "I.mabnn",
                     "V.sa","V.sb","V.sr","T.s",
                     "I.sa0", "I.san",
                     "I.sb0", "I.sbn",
                     "I.sr0", "I.srn",
                     "I.sab00",
                     "I.sar00",
                     "I.sbr00",
                     "I.sabnn",
                     "I.sarnn",
                     "I.sbrnn")
}
names(state) <- state.names

prod.wrapper.a <- function(i,j,n,r = 1, wl = withlike){
    production.fn.a(i,j,n,r,wl)
}
prod.wrapper.b <- function(i,j,n,r = 1, wl = withlike){
    production.fn.b(i,j,n,r,wl)
}
prod.wrapper.r <- function(i,j,n,r = 1, wl = withlike){
    production.fn.r(i,j,n,r,wl)
}

## Run model
dat = data.frame(ode(y = state, times = times,
                     func = wv.BTV.coinfection.reassort, parms = parameters,
                     prod.fn.a=prod.wrapper.a,
                     prod.fn.b=prod.wrapper.b,
                     prod.fn.r=prod.wrapper.r))
save(dat,file="coinfection_reassort_Fu_baseline.RData")
load("coinfection_reassort_Fu_baseline.RData")

## generate plot
tiff("../figures/model_output_fixed_EIP_coinfection_reassort_Fu.tif",res=600,
     compression="lzw",height=600*7,width=600*7)
par(mar=c(2.1,4.1,1.1,2.1),
    oma=c(3,0,3,0))
layout(mat = matrix(c(1,1,2,2), nrow = 2,byrow=TRUE))
plot(dat[,1],(dat[,"V.ma"] + dat[,"V.sa"]),lwd = 2, bty = 'n', 
     xlab = '', ylab = 'Viral load', las = 1,log="y",type='l',
     ylim = range(c(10^fu$log10.mean.titre, (dat[,"V.ma"] + dat[,"V.sa"]))))
lines(dat[,1],(dat[,"V.mb"] + dat[,"V.sb"]),col="gray",lty="dashed",lwd=2)
lines(dat[,1],(dat[,"V.sr"]),col="red",lwd=2)
lines(dat[,1],(rowSums(dat[,grepl("V",names(dat))])),col="green",lwd=1)
points(fu$Time.pi..hours.,10^fu$log10.mean.titre)
legend("bottomright",bty="n",legend=c("Strain A","Strain B","Reassortants", "Total"),
       lty=c("solid","dashed","solid","solid"),col=c("black","gray","red","green"),
       lwd=c(2,2,2,1))
plot(dat[,1],dat[,"V.sr"]/rowSums(dat[,grepl("V",names(dat))]),lwd = 2, bty = 'n', 
      xlab = 'Time (Days)', ylab = 'Proportion reassortant', las = 1,type='l',ylim=c(0,1))
mtext("Time (Days)",side=1,line=3,cex=3/4)
points(samal$time,samal$prop_reassortant,pch=1)
points(combined[gap==0]$time, combined[gap==0]$prop_reassortant,pch=2)
legend("topleft",bty="n",legend=c("Samal, 1987","el Hussein, 1989"),
       pch=c(1,2))
dev.off()

## Fit model
## Combine datasets
combined <- rbind(samal,el.hussein)
combined.runtimes <- combined[,.(runtime=max(time)),by=gap]

## Likelihood function
LL <- function(par) {
    withlike <- par[1]
    quad.a <- par[2]
    quad.b <- par[3]
    quad.c <- par[4]
    prod.wrapper.a <- function(i,j,n,r = 1, wl = withlike){
        production.fn.a(i,j,n,r,wl)
    }
    prod.wrapper.b <- function(i,j,n,r = 1, wl = withlike){
        production.fn.b(i,j,n,r,wl)
    }
    prod.wrapper.r <- function(i,j,n,r = 1, wl = withlike){
        production.fn.r(i,j,n,r,wl)
    }
    loglik <- 0
    for (second.intro in unique(combined.runtimes$gap)){
        runtime <- combined.runtimes[gap==second.intro,runtime]
        print(second.intro)
        if (second.intro > 0) {
            times.1 <- seq(0,second.intro,1)
            times.2 <- seq(second.intro,runtime,1)
            state.1 <- state
            state.1["V.mb"] <- 0
            dat.1 = data.frame(ode(y = state.1, times = times.1,
                                   func = wv.BTV.coinfection.reassort, parms = parameters,
                                   prod.fn.a=prod.wrapper.a,
                                   prod.fn.b=prod.wrapper.b,
                                   prod.fn.r=prod.wrapper.r))
            state.2 <- pmax(as.numeric(dat.1[nrow(dat.1),-1]),0)
            names(state.2) <- names(state)
            state.2["V.mb"] <- state["V.mb"]
            dat.2 = data.frame(ode(y = state.2, times = times.2,
                                   func = wv.BTV.coinfection.reassort, parms = parameters,
                                   prod.fn.a=prod.wrapper.a,
                                   prod.fn.b=prod.wrapper.b,
                                   prod.fn.r=prod.wrapper.r))
            dat <- rbind(dat.1,dat.2)
        } else {
            times <- seq(0,runtime,1)
            dat = data.frame(ode(y = state, times = times,
                                 func = wv.BTV.coinfection.reassort, parms = parameters,
                                 prod.fn.a=prod.wrapper.a,
                                 prod.fn.b=prod.wrapper.b,
                                 prod.fn.r=prod.wrapper.r))
        }
        setDT(dat)
        ## Get probability that BOTH break through (OR a reassortant breaks through)
        p.breakthrough <- 1/(1+exp(-(quad.a * second.intro^2 + quad.b * second.intro + quad.c)))
        
        ## Calculate likelihood
        data <- combined[gap==second.intro, .(time,N_17,N_10,N_reassortant,total)]
        dat[,total.virus:=V.ma+V.mb+V.sa+V.sb+V.sr]
        ## virus.prop <- dat[,.(time,
        ##                      prop.a=(V.ma+V.sa)/total.virus,
        ##                      prop.b=(V.mb+V.sb)/total.virus,
        ##                      prop.r=(V.sr)/total.virus)]
        virus.prop <- dat[,.(time,
                             prop.r=pmax(0,(V.sr)/total.virus))]
        for (ii in 1: nrow(data)) {
            tt <- as.numeric(data[ii,time])
            ## loglik <- loglik + dmultinom(as.numeric(data[ii,.(N_17,N_10,N_reassortant)]),
            ##                              as.numeric(data[ii,.(total)]),
            ##                              prob=as.numeric(virus.prop[time==tt,
            ##                                                         .(prop.a,prop.b,prop.r)]),
            ##                              log=TRUE)
            loglik <- loglik + dbinom(as.numeric(data[ii,.(N_reassortant)]),
                                      as.numeric(data[ii,.(total)]),
                                      prob=p.breakthrough*as.numeric(virus.prop[time==tt,
                                                                                .(prop.r)]),
                                      log=TRUE)
        }
    }
    return(loglik)
}
NLL <- function(par) {
    par[1] <- plogis(par[1])
    par[2] <- -exp(par[2])
    par[3] <- exp(par[3])
    -LL(par)
}

## Try with optim
optim.out <- optim(c(0,0,0,0),NLL)


### Credible and prediction intervals
n.samples <- 100
viral.load <- matrix(NA,n.samples,length(times))
reassortant.load <- viral.load
indices <- sample(nrow(samples),n.samples,replace=TRUE)
for (ii in 1:length(indices)) {
    print(ii)
    index <- indices[ii]
    growth.ratio <- samples[index,1] 
    withlike <- samples[index,2]
    parameters = c(c.m = c.m, c.s = c.s, beta.m = beta.m, beta.s = beta.s, n = n,
                   p.m = p.m, p.s = p.s, epsilon = epsilon, tau = samples[index,3],k = k)
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
    viral.load[ii,] <- rowSums(dat[,grepl("V",names(dat))])
    reassortant.load[ii,] <- dat[,"V.sr"]
}
save(viral.load,reassortant.load,file="coinfection_reassort_Fu_CI.RData")

## generate plot
## tiff("../figures/model_output_fixed_EIP_coinfection_reassort_Fu_CI.tif",res=600,
##      compression="lzw",height=600*7,width=600*7)
plot(dat[,1],apply(reassortant.load/viral.load,2,mean),lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Proportion reassortant', las = 1,ylim=c(0,1),type='l')
lines(dat[,1],apply(reassortant.load/viral.load,2,function(x)quantile(x,0.025)),lty="dotted")
lines(dat[,1],apply(reassortant.load/viral.load,2,function(x)quantile(x,1-0.025)),lty="dotted")
lines(dat[,1],apply(reassortant.load/viral.load,2,function(x)quantile(x,0.25)),lty="dashed",lwd=1.5)
lines(dat[,1],apply(reassortant.load/viral.load,2,function(x)quantile(x,1-0.25)),lty="dashed",lwd=1.5)
mtext("Time (Days)",side=1,line=3,cex=3/4)
points(samal$time,samal$prop_reassortant,pch=1)
points(combined[gap==0]$time, combined[gap==0]$prop_reassortant,pch=2)
legend("topleft",bty="n",legend=c("Samal, 1987","el Hussein, 1989"),
       pch=c(1,2))
## dev.off()



### Different times
for (second.intro in ((1:5)*2-1)*24){
    times.1 <- seq(0,second.intro,1)
    times.2 <- seq(second.intro,20*24,1)
    state.1 <- state
    state.1["V.mb"] <- 0
    dat.1 = data.frame(ode(y = state.1, times = times.1,
                           func = wv.BTV.coinfection.reassort, parms = parameters,
                           prod.fn.a=prod.wrapper.a,
                           prod.fn.b=prod.wrapper.b,
                           prod.fn.r=prod.wrapper.r))
    state.2 <- as.numeric(dat.1[nrow(dat.1),-1])
    names(state.2) <- names(state)
    state.2["V.mb"] <- initial.titre/2
    dat.2 = data.frame(ode(y = state.2, times = times.2,
                           func = wv.BTV.coinfection.reassort, parms = parameters,
                           prod.fn.a=prod.wrapper.a,
                           prod.fn.b=prod.wrapper.b,
                           prod.fn.r=prod.wrapper.r))
    dat <- rbind(dat.1,dat.2)
    tiff(paste0("../figures/model_output_fixed_EIP_coinfection_reassort_sep",second.intro,
                "_fitted_Fu.tif"),
         res=600,compression="lzw",height=600*7,width=600*7)
    par(mar=c(2.1,4.1,1.1,2.1),
        oma=c(3,0,3,0))
    layout(mat = matrix(c(1,1,2,2), nrow = 2,byrow=TRUE))
    plot(dat[,1],(dat[,"V.ma"] + dat[,"V.sa"]),lwd = 2, bty = 'n', 
         xlab = '', ylab = 'Viral load', las = 1,log="y",type='l',
         ylim = range(c(10^fu$log10.mean.titre, (dat[,"V.ma"] + dat[,"V.sa"]))))
    lines(dat[,1],(dat[,"V.mb"] + dat[,"V.sb"]),col="gray",lty="dashed",lwd=2)
    lines(dat[,1],(dat[,"V.sr"]),col="red",lwd=2)
    lines(dat[,1],(rowSums(dat[,grepl("V",names(dat))])),col="green",lwd=1)
    points(fu$Time.pi..hours.,10^fu$log10.mean.titre)
    legend("bottomright",bty="n",legend=c("Strain A","Strain B","Reassortants", "Total"),
           lty=c("solid","dashed","solid","solid"),col=c("black","gray","red","green"),
           lwd=c(2,2,2,1))
    plot(dat[,1],dat[,"V.sr"]/rowSums(dat[,grepl("V",names(dat))]),lwd = 2, bty = 'n', 
         xlab = 'Time (Days)', ylab = 'Proportion reassortant', las = 1,type='l',ylim=c(0,1))
    mtext("Time (Days)",side=1,line=3,cex=3/4)
    points(combined[gap==second.intro]$time,
           combined[gap==second.intro]$prop_reassortant,pch=2)
    legend("topleft",bty="n",legend=c("el Hussein, 1989"),
           pch=c(2),col=c("black"))
    dev.off()
}


### Explore parameter space - NEEDS UPDATING
growth.ratios <- 10^seq(0,10)
withlikes <- seq(0,1,0.1)
param.grid <- expand.grid(growth.ratios,withlikes)
output.grid <- vector(mode="numeric",length=nrow(param.grid))

for (iii in 1:nrow(param.grid)) {
    print(iii)
    growth.ratio.temp <- param.grid[iii,1]
    withlike.temp <- param.grid[iii,2]
    prod.wrapper.a <- function(i,j,n,r = growth.ratio.temp, wl = withlike.temp){
        production.fn.a(i,j,n,r,wl)
    }
    prod.wrapper.b <- function(i,j,n,r = growth.ratio.temp, wl = withlike.temp){
        production.fn.b(i,j,n,r,wl)
    }
    prod.wrapper.r <- function(i,j,n,r = growth.ratio.temp, wl = withlike.temp){
        production.fn.r(i,j,n,r,wl)
    }
    dat = data.frame(ode(y = state, times = times,
                         func = wv.BTV.coinfection.reassort, parms = parameters,
                         prod.fn.a=prod.wrapper.a,
                         prod.fn.b=prod.wrapper.b,
                         prod.fn.r=prod.wrapper.r))
    output.grid[iii] <- dat[nrow(dat),"V.sr"]/sum(dat[nrow(dat),grepl("V",names(dat))])
}

## save(output.grid, file="../output/output_grid.RData")
load("../output/output_grid.RData")

## Contour plot
xtick <- pretty(growth.ratios)
x <- log(growth.ratios,base=10)
ytick <- pretty(withlikes)
tiff(paste0("../figures/contour_plot_growth_ratio_withlike.tif"),
     res=600,compression="lzw",height=600*7,width=600*7)
filled.contour(x, withlikes,
               matrix(output.grid,nrow=length(growth.ratios)),
               color.palette=magma,
               plot.axes = {axis(1,at=pretty(x),label=10^pretty(x));axis(2)},
               xlab="Growth ratio",
               ylab="With-like proportion")
points(log(growth.ratio,base=10), withlike,col="red",pch=19)
text(1.1*log(growth.ratio,base=10), withlike,"Best-fitting\nparameter combination",adj=0,col="red")
dev.off()


### barriers - basic
##rerun
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

## barrier stuff
barrier.prob <- 0.3
secondary.barrier.prob <- 12/30
infection.prob <- barrier.prob*(2 - barrier.prob)
load("coinfection_reassort_Fu_baseline.RData")
reassort.prop <- dat[,"V.sr"]/rowSums(dat[,grepl("V",names(dat))])
n.plaques <- 20
pi.zeroes <- 2*(1-barrier.prob)/(2-barrier.prob)
## generate plot
tiff("../figures/model_output_fixed_EIP_coinfection_reassort_barrier_basic_Fu.tif",res=600,
     compression="lzw",height=600*7,width=600*7)
par(mar=c(2.1,4.1,1.1,2.1),
    oma=c(3,0,3,0))
## layout(mat = matrix(c(1,1,2,2), nrow = 2,byrow=TRUE))
## plot(dat[,1],0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]),lwd = 2, bty = 'n', 
##      xlab = '', ylab = 'Viral load', las = 1,log="y",type='l',
##      ylim = range(c(10^mellor$titre, 0.5*initial.titre*(dat[,"V.ma"] + dat[,"V.sa"]))))
## lines(dat[,1],0.5*initial.titre*(dat[,"V.mb"] + dat[,"V.sb"]),col="gray",lty="dashed",lwd=2)
## lines(dat[,1],0.5*initial.titre*(dat[,"V.sr"]),col="red",lwd=2)
## lines(dat[,1],0.5*initial.titre*(rowSums(dat[,grepl("V",names(dat))])),col="green",lwd=1)
## points(mellor$day,10^mellor$titre)
## legend("topleft",bty="n",legend=c("Strain A","Strain B","Reassortants", "Total"),
##        lty=c("solid","dashed","solid","solid"),col=c("black","gray","red","green"),
##        lwd=c(2,2,2,1))
plot(dat[,1],reassort.prop,
     lwd = 2, bty = 'n', xaxs="i",yaxs="i",
     xlab = 'Time (Days)', ylab = 'Proportion reassortant', las = 1,type='l',ylim=c(0,1))
## polygon(c(dat[,1],rev(dat[,1])),
##         c(reassort.upper,rev(reassort.lower)),
##         col=adjustcolor("red",alpha.f=0.25),border=FALSE)
lines(dat[,1],qzib(0.975,n.plaques,prob=reassort.prop,pi=pi.zeroes)/n.plaques,
      col="red",lty="dotted")
lines(dat[,1],qzib(0.025,n.plaques,prob=reassort.prop,pi=pi.zeroes)/n.plaques,
      col="red",lty="dotted")
lines(dat[,1],qzib(0.95,n.plaques,prob=reassort.prop,pi=pi.zeroes)/n.plaques,
      col="red",lty="dashed")
lines(dat[,1],qzib(0.05,n.plaques,prob=reassort.prop,pi=pi.zeroes)/n.plaques,
      col="red",lty="dashed")
lines(dat[,1],qzib(0.75,n.plaques,prob=reassort.prop,pi=pi.zeroes)/n.plaques,
      col="red")
lines(dat[,1],qzib(0.25,n.plaques,prob=reassort.prop,pi=pi.zeroes)/n.plaques,
      col="red")
mtext("Time (Days)",side=1,line=3,cex=3/4)
points(samal$time,samal$prop_reassortant,pch=1)
points(combined[gap==0]$time, combined[gap==0]$prop_reassortant,pch=2)
legend("topleft",bty="n",legend=c("Samal, 1987","el Hussein, 1989"),
       pch=c(1,2))
dev.off()





### BONEYARD
### With temperature dependence - NEEDS WORK!!!!!!!
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
points(combined[gap==0]$time, combined[gap==0]$prop_reassortant,pch=2)
points(0.95*max(dat[,1]),
       mean(c(samal$prop_reassortant,combined[gap==0]$prop_reassortant)),
       col="red",pch=23)
legend("topleft",bty="n",legend=c("Samal, 1987","el Hussein, 1989", "Mean"),
       pch=c(1,2,23),col=c("black","black","red"))
dev.off()

## get quantiles
CI.barrier.fn <- function(upper=0.975,lower=0.025,n.plaques=20) {
    if (barrier.prob^2/infection.prob > upper) {
        reassort.upper <- qbinom((upper + barrier.prob^2/infection.prob - 1)/(barrier.prob^2/infection.prob),
                                 n.plaques,reassort.prop)/n.plaques
        reassort.lower <- qbinom((lower + barrier.prob^2/infection.prob - 1)/(barrier.prob^2/infection.prob),
                                 n.plaques,reassort.prop)/n.plaques
    } else if (barrier.prob^2/infection.prob > lower) {
        reassort.upper <- qbinom((upper + barrier.prob^2/infection.prob - 1)/(barrier.prob^2/infection.prob),
                                 n.plaques,reassort.prop)/n.plaques
        reassort.lower <- rep(0,nrow(dat))
    } else {
        reassort.upper <- rep(0,nrow(dat))
        reassort.lower <- rep(0,nrow(dat))
    }
    return(list(upper=reassort.upper,lower=reassort.lower))
}
