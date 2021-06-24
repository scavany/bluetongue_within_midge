## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspacer
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(data.table,
       viridis,
       fitR,
       deSolve,
       grDevices
       ##,BayesianTools,psych,bbmle
       )

## load within-vector model of BTV infection
source("./within_midge_barrier_based_fn.R")

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
initial.titre <- 10^fu$log10.mean.titre[1] 

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)

## Set up timing vectors
## times = seq(from = 0, to = 14*24, by = 1)
## data.hours.it <- which(times %in% fu.intrathoracic.hours)
## model.hours.it <- which(fu.intrathoracic.hours %in% times)
## titre.data.it <- floor(fu.intrathoracic$titre[model.hours.it]+0.5)
## data.hours <- which(times %in% fu$Time.pi..hours.)
## model.hours <- which(fu$Time.pi..hours. %in% times)
## titre.data <- floor(10^fu$log10.mean.titre[model.hours] + 0.5)

## Barrier probabilities
barrier.probs <- c(p.mib=1/3,p.meb=1/2,p.db=1/5,
                   p.sgib=1,p.sgeb=1,p.totb=0)

## Oral infection - parms and initial conditions 
parms <- c(c.b=3, beta.b=1e-2,
           c.m=3, d.m = 0.1, p.m=1e3, beta.m=1e-2,
           c.l=0.3, d.l = 0.1, p.l=1e2, beta.l=1e-2,
           c.d=0.3, d.d = 0.1, p.d=1e2, beta.d=1e-2,
           c.s=0.3, p.s=1e2, beta.s=1e-2)

state <- c(V.b = floor(initial.titre*log(2) + 0.5), V.m=0, V.l=0, V.d=0, V.s=0,
           N.m=1e3, N.l=1e3, N.d=1e3, N.s=1e2,
           T.m=0, T.l=0, T.d=0, T.s=0,
           I.m=0, I.l=0, I.d=0, I.s=0)

times <- seq(0,10,1/24)

## Run stochastic model
n.sims <- 1000
trajs <- array(NA,dim=c(n.sims,dim(wv.BTV.barrier.stoch(times,state,parms,barrier.probs))),
               dimnames=list(1:n.sims,
                             rownames(wv.BTV.barrier.stoch(times,state,parms,barrier.probs)),
                             colnames(wv.BTV.barrier.stoch(times,state,parms,barrier.probs))))

## trajs <- list()
for (iii in 1:n.sims) {
    ## print(iii)
    ## trajs[[iii]] <- wv.BTV.barrier.stoch(times,state,parms,barrier.probs)
    trajs[iii,,] <- as.matrix(wv.BTV.barrier.stoch(times,state,parms,barrier.probs))
}

## Plot
pdf("../figures/barrier_model_plot.pdf",
    width=14,height=7,pointsize=10)
par(mfrow=c(1,2))
plot(-1,0.1,xlim=c(0,10),ylim=c(1,1e6),
     xlab="Time (days)",ylab="Virions",log="y",bty="n",xaxs="i")
cols.alpha <- adjustcolor(viridis(5),0.4)
cols <- viridis(5)
for (iii in 1:n.sims) {
    ## print(iii)
    lines(trajs[iii,,"time"],trajs[iii,,"V.b"],col=cols.alpha[1])
    lines(trajs[iii,,"time"],trajs[iii,,"V.m"],col=cols.alpha[2])
    lines(trajs[iii,,"time"],trajs[iii,,"V.l"],col=cols.alpha[3])
    lines(trajs[iii,,"time"],trajs[iii,,"V.d"],col=cols.alpha[4])
    lines(trajs[iii,,"time"],trajs[iii,,"V.s"],col=cols.alpha[5])
}
lines(trajs[1,,"time"],apply(trajs[,,"V.b"],c(2),mean),
      col=cols[1],lwd=3)
lines(trajs[1,,"time"],apply(trajs[,,"V.m"],c(2),mean),
      col=cols[2],lwd=3)
lines(trajs[1,,"time"],apply(trajs[,,"V.l"],c(2),mean),
      col=cols[3],lwd=3)
lines(trajs[1,,"time"],apply(trajs[,,"V.d"],c(2),mean),
      col=cols[4],lwd=3)
lines(trajs[1,,"time"],apply(trajs[,,"V.s"],c(2),mean),
      col=cols[5],lwd=3)
abline(h=10^0.75,lty="dotted",lwd=2)
abline(h=10^3,lty="dashed",lwd=2)
legend("topleft",bty="n",legend=c("Blood meal","Midgut","Fat-body cells","Ganglia","Salivary glands"),
       lty=1,col=viridis(5),lwd=2)
plot(-1,0.1,xlim=c(0,10),ylim=c(1,1e6),
     xlab="Time (days)",ylab="Virions",log="y",bty="n",xaxs="i")
trajs.tot <- apply(trajs[,,c("V.b","V.m","V.l","V.d","V.s")],c(1,2),sum)
trajs.pos <- trajs.tot
trajs.pos[which(trajs.tot < 10^0.75)] <- NA
for (iii in 1:n.sims) {
    ##print(iii)
    lines(trajs[iii,,"time"],trajs.tot[iii,])
}
lines(trajs[1,,"time"],
      ## rowSums(apply(trajs[,,c("V.b","V.m","V.l","V.d","V.s")],c(2,3),mean)),
      apply(trajs.pos,2,function(x)mean(x,na.rm=T)),
      col="red",lwd=3)
abline(h=10^0.75,lty="dotted",lwd=2)
abline(h=10^3,lty="dashed",lwd=2)
points(fu$Time.pi..hours/24,10^fu$log10.mean.titre,col="green",pch=19)
legend("topleft",bty="n",legend=c("Individual midges","Mean of positive midges","Data (Fu, 1999)"),
       lty=c(1,1,NA),col=c("black","red","green"),lwd=c(1,3,NA),pch=c(NA,NA,19))
dev.off()

