## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspacer
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(data.table,
       viridis,
       deSolve,
       grDevices,
       BayesianTools,
       mixtools,
       adaptivetau,
       RColorBrewer,
       imager)

## load within-vector model of BTV infection
source("./within_midge_barrier_based_fn.R")

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
initial.titre <- 10^fu$log10.mean.titre[1]
equilibrium.titre <- mean(tail(10^fu$log10.mean.titre)) ##Should this be the mean of the log10??
lod <- 10^0.75

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)
initial.titre.it <- fu.intrathoracic$titre[1]
equilibrium.titre.it <- mean(tail(fu.intrathoracic$titre)) ##Should this be the mean of the log10??

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
pdf("../figures/logistic_prevalence.pdf")
plot(times/24,prop.positive.fn(times,logistic.optim.out$par[1],
                               logistic.optim.out$par[2],
                               logistic.optim.out$par[3]),type='l',
     ylim=c(0,1.05),xlim=c(0,15),
     ylab="Proportion BTV positive",xlab="Days post-infection",las=1,xaxs="i",yaxs="i",
     bty="n",lwd=3)
points(fu$Time.pi..hours./24,fu$Detection.rate)
legend("topright",pch=c(1,NA),lty=c(NA,1),legend=c("Data from Fu et al.","Fitted estimate"),
       bty="n",lwd=c(NA,3))
dev.off()
prop.positive.vec <- prop.positive.fn(times,logistic.optim.out$par[1],
                                      logistic.optim.out$par[2],
                                      logistic.optim.out$par[3])
prop.disseminated <- min(prop.positive.vec)

## Get expression for cb and cm
c.m <- log(initial.titre/lod) / logistic.optim.out$par[3]

## Barrier probabilities
p.mib=1/3;p.db=0.12
barrier.probs <- c(p.mib=p.mib,p.db=p.db)

## Oral infection - parms and initial conditions
parms <- c(c.b=c.m, beta.b=1e-5,
           c.m=c.m, beta.m=1e-5, d.m=0.1, p.m=100, T0.m=100, mu.m=1, mui.m=0.1,epsilon.m=0,
           c.s=1, beta.s=1e-5, d.s=0.1, p.s=1000, T0.s=1000, mu.s=1, mui.s=0.1,epsilon.s=0)

state <- c(V.b = floor(initial.titre*log(2) + 0.5),
           V.m=0, T.m=0, E.m = 0, I.m=0, 
           V.s=0, T.s=0, E.s = 0, I.s=0)
state <- c(V.b = floor(initial.titre*log(2) + 0.5),
           V.m=0, T.m=0, E.m = 0, E.mm = 0, I.m=0, I.mm=0, 
           V.s=0, T.s=0, E.s = 0, E.ss = 0, I.s=0, I.ss=0
           )

### Load mcmc results
##load(file="./barrier_mcmc_out.RData",verbose=T) # fit to all data
##load(file="./barrier_mcmc_out_oral_only.RData",verbose=T) # fit only to oral data
##load(file="./barrier_mcmc_out_backflow.RData",verbose=T) # allow backflow
##load(file="./barrier_mcmc_out_backflow_oral_only.RData",verbose=T) # fit to oral only, allow backflow
##load(file="./barrier_mcmc_out_extra_constraints.RData",verbose=T) # extra constraints
load(file="./coinfection_mcmc_out_joint_fitting_10000.RData",verbose=T) # joint fit

pdf("../figures/mcmc_out.pdf")
plot(mcmc.out,start=mcmc.out$settings$iterations/9)
dev.off()

##pdf("../figures/corr_plot.pdf")
##correlationPlot(mcmc.out,start=mcmc.out$settings$iterations / 9)
##dev.off()
sample.out <- as.data.frame(getSample(mcmc.out,start=mcmc.out$settings$iterations / 9))
MAP.out <- MAP(mcmc.out,start=mcmc.out$settings$iterations / 9)$parametersMAP

### Simple parameter exploration
## contour of probability of zero infections, mechanism: bloodmeal depletion before infection,
## but with an MIB
beta.T.m <- 10^seq(-4,1,0.1)
n.virions <- 10^seq(0,log(10*initial.titre*log(2),10),0.1)
prob.0 <- expand.grid(n.virions=n.virions,beta.T.m=beta.T.m)
prob.0$prob.0 <- p.mib * (1 - (c.m / (c.m + prob.0$beta.T.m)) ^ prob.0$n.virions)
zmat <- matrix(prob.0$prob.0,nrow=length(n.virions),ncol=length(beta.T.m))
mode.out <- MAP.out["beta.mT.m"]#mode.out <- MAP.out["beta.m"]*MAP.out["T0.m"]
q.5 <- quantile(sample.out$beta.mT.m,c(0.25,0.75))
q.05 <- quantile(sample.out$beta.mT.m,c(0.025,0.975))
## plot
tiff("../figures/midgut_infection_prob_fixed_barrier.tif",
     width=1200,height=900,res=300,compression="lzw",pointsize=10)
par(mar=c(4.1,4.1,2,0.1))
filled.contour(log(n.virions,10),log(beta.T.m,10),zmat,
               zlim=c(0,1),
               plot.title = title(#"Probability that at least one midgut cell is infected",
                                  xlab="Initial number of virions",
                                  ylab=expression(paste(beta,T[m]," (infected cells / hr / virion)"))),
               axes=FALSE,
               plot.axes={axis(1,at=pretty(range(log(n.virions,10)),6),
                               labels=round(10^pretty(range(log(n.virions,10)),6),2));
                               axis(2,at=pretty(range(log(beta.T.m,10)),5),
                                    labels=round(10^pretty(range(log(beta.T.m,10)),5),5));
                               abline(v=log(initial.titre*log(2),10),lty=2,lwd=2);
                               points(log(initial.titre*log(2),10),log(mode.out,10),
                                      pch=20,cex=2);
                               lines(c(0.95,1.05)*log(initial.titre*log(2),10),
                                     log(rep(q.5[1],2),10));
                               lines(c(0.95,1.05)*log(initial.titre*log(2),10),
                                     log(rep(q.5[2],2),10));
                               lines(c(0.95,1.05)*log(initial.titre*log(2),10),
                                     log(rep(q.05[1],2),10),lty=3);
                               lines(c(0.95,1.05)*log(initial.titre*log(2),10),
                                     log(rep(q.05[2],2),10),lty=3);
               },
               key.axes=axis(4,at=pretty(c(zmat,1))),
               key.title={par(cex.main=0.9);title("Pr(midgut\ninfection)")},
               color=colorRampPalette(brewer.pal(9,"Reds")))
dev.off()

## as above, but vary p.mib
beta.T.m <- 10^seq(-4,1,0.1)
n.virions <- log(2) * initial.titre
p.mib.vec <- seq(0,1,0.01)
prob.0 <- expand.grid(p.mib=p.mib.vec,beta.T.m=beta.T.m)
prob.0$prob.0 <- prob.0$p.mib * (1 - (c.m / (c.m + prob.0$beta.T.m)) ^ n.virions)
zmat <- matrix(prob.0$prob.0,nrow=length(p.mib.vec),ncol=length(beta.T.m))
## plot
tiff("../figures/midgut_infection_prob_fixed_bloodmeal.tif",
     width=1200,height=900,res=300,compression="lzw",pointsize=10)
par(mar=c(4.1,4.1,2,0.1))
filled.contour(p.mib.vec,log(beta.T.m,10),zmat,
               plot.title = title(#"Probability that at least one midgut cell is infected",
                                  xlab="Probability not an MIB",
                                  ylab=expression(paste(beta,T[m]," (infected cells / hr / virion)"))),
               axes=FALSE,
               plot.axes={axis(1,at=pretty(range(p.mib.vec),6),
                               labels=round(pretty(range(p.mib.vec),6),2));
                               axis(2,at=pretty(range(log(beta.T.m,10)),5),
                                    labels=round(10^pretty(range(log(beta.T.m,10)),5),5));
                               abline(v=p.mib,lty=2,lwd=2);
                               points(p.mib,log(mode.out,10),
                                      pch=20,cex=2);
                               lines(c(0.9,1.1)*p.mib,
                                     log(rep(q.5[1],2),10));
                               lines(c(0.9,1.1)*p.mib,
                                     log(rep(q.5[2],2),10));
                               lines(c(0.9,1.1)*p.mib,
                                     log(rep(q.05[1],2),10),lty=3);
                               lines(c(0.9,1.1)*p.mib,
                                     log(rep(q.05[2],2),10),lty=3);
               },
               key.axes=axis(4,at=pretty(zmat)),
               key.title={par(cex.main=0.9);title("Pr(midgut\ninfection)")},
               color=colorRampPalette(brewer.pal(9,"Reds")))
dev.off()

## as above, but vary p.mib and bloodmeal, fixing parms to MAP
beta.T.m <- as.numeric(mode.out)
n.virions <- 10^seq(0,log(10*initial.titre*log(2),10),0.1)
prob.0 <- expand.grid(p.mib=p.mib.vec,n.virions=n.virions)
prob.0$prob.0 <- prob.0$p.mib * (1 - (c.m / (c.m + beta.T.m)) ^ prob.0$n.virions)
zmat <- matrix(prob.0$prob.0,nrow=length(p.mib.vec),ncol=length(n.virions))
## plot
tiff("../figures/midgut_infection_prob_fixed_parms.tif",
     width=1200,height=900,res=300,compression="lzw",pointsize=10)
par(mar=c(4.1,4.1,2,0.1))
filled.contour(p.mib.vec,log(n.virions,10),zmat,
               plot.title = title(#"Probability that at least one midgut cell is infected",
                                  xlab="Probability not an MIB",
                                  ylab="Initial number of virions"),
               axes=FALSE,
               plot.axes={axis(1,at=pretty(range(p.mib.vec),6),
                               labels=round(pretty(range(p.mib.vec),6),2));
                               axis(2,at=pretty(range(log(n.virions,10)),5),
                                    labels=round(10^pretty(range(log(n.virions,10)),5),5));
                               abline(v=p.mib,lty=2,lwd=2);
                               abline(h=log(log(2) * initial.titre,10),lty=2,lwd=2)
               },
               key.axes=axis(4,at=pretty(zmat)),
               key.title={par(cex.main=0.9);title("Pr(midgut\ninfection)")},
               color=colorRampPalette(brewer.pal(9,"Reds")))
dev.off()

pdf("../figures/combined_barrier_plot.pdf",width=7,height=3.5)
par(mfrow=c(1,2))
par(mar = c(0,0,0,0))
im <- load.image("../figures/midgut_infection_prob_fixed_bloodmeal.tif")
plot(im,axes=FALSE)
mtext(side = 3, line = -1.2, adj = 0.04, 'A', font = 2)
im <- load.image("../figures/midgut_infection_prob_fixed_barrier.tif")
plot(im,axes=FALSE)
mtext(side = 3, line = -1.2, adj = 0.04, 'B', font = 2)
##par(mar = c(0,0.1,1.3,0.1))
dev.off()

## Get parms
parms <- MAP.out
parms["beta.m"] <- parms[["beta.mT.m"]]/parms[["T0.m"]]
parms["c.m"] <- c.m;parms["c.b"] <- c.m;parms["beta.b"] <- parms["beta.m"];parms["mui.m"] <- parms["mu.m"];parms["mui.s"] <- parms["mu.s"]
V.out <- calc.V(state,parms,times)
V.adj <- (prop.positive.vec - prop.disseminated) * V.out$V.tot + prop.disseminated * V.out$V.tot.pos
V.tot.it <- calc.V.it(state,parms,times)$V.tot.it


pdf("../figures/viral_dynamics_det.pdf")
cols <- plasma(6)[c(1,4)]
par(mfrow=c(2,1),mar=c(5.1,4.1,1.1,2.1))
plot(times, V.adj / log(2),type='l', lwd=3, las=1,xaxs="i",yaxs="i",log="y",bty="n",col=cols[1],
     xlab="", ylab=expression(Viral~load~(TCID[50])),
     ylim=c(min(c(V.adj / log(2),titre.data,0.8*10^0.75)),max(c(V.out$V.tot.3 / log(2),titre.data))))
lines(times,V.out$V.tot.1 / log(2),col=cols[2],lwd=2,lty=2)
lines(times,V.out$V.tot.2 / log(2),col="gray",lwd=2,lty=2)
lines(times,V.out$V.tot.3 / log(2),col="black",lwd=2,lty=2)
legend("bottomright",
       legend=c("All positive midges (model)","Midges with an MIB (model)","Midges without an MIB (model)",
                "Midges with no barriers (model)","Average of positive midges (data)"),
       lwd=c(3,2,2,2,1),pch=c(rep(NA,4),1),
       lty=c(1,2,2,2,NA),col=c(cols,"gray","black","black"),bty="n")
points(times[data.hours],titre.data)
abline(h=10^3, lty=2)
abline(h=10^0.75, lty=3)
mtext(side = 3,'A',adj=0)
plot(times, V.tot.it / log(2),type='l', lwd=3, las=1,xaxs="i",yaxs="i",log="y",bty="n",
     xlab="Time (hours)", ylab=expression(Viral~load~(TCID[50])),col=cols[1],
     ylim=c(min(c(V.tot.it / log(2),titre.data.it)),max(c(V.tot.it / log(2),titre.data.it))))
points(times[data.hours.it],titre.data.it)
mtext(side = 3,'B',adj=0)
dev.off()

## Plot epsilon against mean EIP, holding other parameters fixed
## temperature function based on Wittmann (default BTV-10)
temperature <- function(EIP,a=0.0069,b=-0.0636) {
    (1 / EIP - b) / a
}
EIP <- function(temperature,a=0.0069,b=-0.0636) {
    1 / (a * temperature + b)
}

## Get epsilon - EIP relationship
times <- seq(0,24*EIP(15),1)
prop.positive.vec <- prop.positive.fn(times,logistic.optim.out$par[1],
                                      logistic.optim.out$par[2],
                                      logistic.optim.out$par[3])
epsilon.ms <- exp(seq(log(parms["epsilon.m"]/sqrt(10)),
                      log(parms["epsilon.m"]*sqrt(10)),
                      length.out=201))
epsilon.ss <- exp(seq(log(parms["epsilon.s"]/sqrt(10)),
                      log(parms["epsilon.s"]*sqrt(10)),
                      length.out=length(epsilon.ms)))
parms.temp <- parms
eip.vec <- rep(NA,length(epsilon.ms))
for (ii in 1:length(epsilon.ms)) {
    parms.temp["epsilon.m"] <- epsilon.ms[ii]
    parms.temp["epsilon.s"] <- epsilon.ss[ii]
    V.out <- calc.V(state,parms.temp,times)
    V.adj <- (prop.positive.vec - prop.disseminated) * V.out$V.tot + prop.disseminated * V.out$V.tot.pos
    min.titre.time <- times[which.min(V.adj)]
    if (length(times[which(V.adj > 10^2.5 * log(2) & times > min.titre.time)])) {
        eip.vec[ii] <- min(times[which(V.adj > 10^2.5 * log(2) & times > min.titre.time)])
    }
}

## Plot it
png("../figures/epsilon_eip_relationship.png")
par(mar=c(5.1,4.1,2.1,4.1))
cols <- brewer.pal(n=3,name="Reds")
plot(epsilon.ms,eip.vec,log="xy",type='l',xlab="",
     ylab="Extrinsic incubation period (hours)",
     bty="n",las=1,#ylim=c(1,EIP(15))*24,
     yaxs="i",lwd=2.5,xaxt="n")
## abline(h=max(times))
abline(v=epsilon.ms[length(eip.vec)/2 + 0.5],lty="dashed")
abline(h=eip.vec[length(eip.vec)/2 + 0.5],lty="dashed")
axis(1)
axis(1,at=approxfun(epsilon.ss,epsilon.ms)(axTicks(1)),
     labels=FALSE,#axTicks(1),
     col.ticks=cols[3],col.axis=cols[3])
mtext(expression(epsilon[m]*", "),side=1,line=3)
mtext(expression("\t"*epsilon[s]),side=1,line=3,col=cols[3])
mtext("\t\t\t (hours)",side=1,line=3)
axis(4,at=axTicks(2),las=1,
     labels=ifelse(is.infinite(temperature(axTicks(2))),expression(infinity),
                   round(temperature(axTicks(2) / 24))))
mtext("Temperature (C) for BTV-10",4,3)
dev.off()

