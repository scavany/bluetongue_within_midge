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
       grDevices,
       BayesianTools,
       mixtools,
       RColorBrewer#,psych,bbmle
       )

## load within-vector model of BTV infection
source("./within_midge_barrier_based_fn.R")

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
initial.titre <- 10^fu$log10.mean.titre[1]
equilibrium.titre <- mean(tail(10^fu$log10.mean.titre)) ##Should this be the mean of the log10??

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
## plot(times,prop.positive.fn(times,logistic.optim.out$par[1],
##                             logistic.optim.out$par[2],
##                             logistic.optim.out$par[3]),type='l',
##      ylim=c(0,1),ylab="Proportion positive")
## points(fu$Time.pi..hours.,fu$Detection.rate)
prop.positive.vec <- prop.positive.fn(times,logistic.optim.out$par[1],
                                      logistic.optim.out$par[2],
                                      logistic.optim.out$par[3])
prop.disseminated <- min(prop.positive.vec)

## Barrier probabilities
p.mib=1/3;p.db=0.12
barrier.probs <- c(p.mib=p.mib,p.db=p.db)

## Oral infection - parms and initial conditions
parms <- c(c.b=1, beta.b=1e-5,
           c.m=1, beta.m=1e-5, d.m=0.1, p.m=100, T0.m=100, mu.m=1, mui.m=0.1,epsilon.m=0,
           c.s=1, beta.s=1e-5, d.s=0.1, p.s=1000, T0.s=1000, mu.s=1, mui.s=0.1,epsilon.s=0)

state <- c(V.b = floor(initial.titre*log(2) + 0.5),
           V.m=0, T.m=0, E.m = 0, I.m=0, 
           V.s=0, T.s=0, E.s = 0, I.s=0)

### Load mcmc results
load(file="./barrier_mcmc_out_only_oral.RData",verbose=T)

pdf("../figures/mcmc_out.pdf")
plot(mcmc.out,start=1.5e6)
dev.off()

pdf("../figures/corr_plot.pdf")
correlationPlot(mcmc.out,start=1.5e6)
dev.off()

sample.out <- as.data.frame(getSample(mcmc.out,start=1.5e6))
MAP.out <- MAP(mcmc.out,start=1.5e6)$parametersMAP

### Simple parameter exploration
## contour of probability of zero infections, mechanism: bloodmeal depletion before infection
## Make this max out at a third, and instead show parameter combinations which certain dose response would give.
## Also show the dose-response curve? For modal parameter set?
c <- 10^seq(-2,1,0.1)
beta.T.m <- 10^seq(-5,-1,0.1)
n.virions <- initial.titre*log(2)
prob.0 <- expand.grid(c=c,beta.T.m=beta.T.m)
prob.0$prob.0 <- (prob.0$c / (prob.0$c + prob.0$beta.T.m)) ^ n.virions
zmat <- matrix(prob.0$prob.0,nrow=length(c),ncol=length(beta.T.m))
## c.m vs beta.m*T.m to give probability of a third
prob.per.virion <- (1/3)^(1/n.virions)
beta.T.m.baseline <- c * (1 - prob.per.virion) / prob.per.virion
temp.df <- data.frame(sample.out$c.m,sample.out$beta.m*sample.out$T0.m)
median.out <- sapply(temp.df,median)
mean.out <- colMeans(temp.df)
mode.out <- c(MAP.out["c.m"],MAP.out["beta.m"]*MAP.out["T0.m"])
cov.out <- cov(temp.df)
ell.25 <- ellipse(mean.out,cov.out,lwd=2,alpha=0.25,draw=FALSE)
ell.05 <- ellipse(mean.out,cov.out,lwd=2,alpha=0.05,draw=FALSE)
## plot
filled.contour(log(c,10),log(beta.T.m,10),zmat,
               plot.title = title("Probability that no midgut cells are infected",
                                  xlab=expression(c[m]),ylab=expression(paste(beta,T[m]))),
               axes=FALSE,
               plot.axes={axis(1,at=pretty(range(log(c,10)),6),
                               labels=round(10^pretty(range(log(c,10)),6),2));
                               axis(2,at=pretty(range(log(beta.T.m,10)),6),
                                    labels=round(10^pretty(range(log(beta.T.m,10)),6),4));
                               lines(log(c,10),log(beta.T.m.baseline,10),lwd=3);
                               lines(log(ell.25,10),lwd=2);
                               lines(log(ell.05,10),lwd=2,lty="dashed");
                               ## points(median.out[1],median.out[2],pch=1)
                               ## points(mean.out[1],mean.out[2],pch=20)
                               points(log(mode.out[1],10),log(mode.out[2],10),pch=20)
                               text((min(log(c,10)) + 3 * max(log(c,10)))/4,
                                    (min(log(beta.T.m.baseline,10))  + 3 * max(log(beta.T.m.baseline,10)))/4,
                                    labels="P = 1/3", font=2,
                                    srt=atan(diff(range(log(beta.T.m.baseline,10)))/diff(range(log(c,10))))*180/pi,
                                    pos=3)},
               key.axes=axis(4,at=pretty(zmat)),
               color=colorRampPalette(brewer.pal(9,"Reds")))

## Plot the R0 as a function of c and p/mu.i, assuming that 1/3 of midgut infections are established
p.mu.i <- 10^seq(0,4,0.1)
R0 <- expand.grid(c=c,p.mu.i=p.mu.i)
R0$beta.T.m <- beta.T.m.baseline[match(R0$c,c)]
R0$R0 <- R0$beta.T.m / (R0$c + R0$beta.T.m) * R0$p.mu.i
zmat <- matrix(R0$R0,nrow=length(c),ncol=length(p.mu.i))
temp.df <- data.frame(sample.out$c.m,sample.out$p.m*sample.out$mu.m)
median.out <- sapply(temp.df,median)
mean.out <- colMeans(temp.df)
mode.out <- c(MAP.out["c.m"],MAP.out["p.m"]*MAP.out["mu.m"])
cov.out <- cov(temp.df)
ell.25 <- ellipse(mean.out,cov.out,lwd=2,alpha=0.25,draw=FALSE)
ell.05 <- ellipse(mean.out,cov.out,lwd=2,alpha=0.05,draw=FALSE)
## plot
filled.contour(log(c,10),log(p.mu.i,10),zmat,
               plot.title = title("Midgut R0, assuming 1/3 of infections establish in midgut",
                                  xlab=expression(c[m]),ylab=expression(paste(p,"/",mu[i]))),
               axes=FALSE,
               plot.axes={axis(1,at=pretty(range(log(c,10)),6),
                               labels=round(10^pretty(range(log(c,10)),6),2));
                               axis(2,at=pretty(range(log(p.mu.i,10)),6),
                                    labels=round(10^pretty(range(log(p.mu.i,10)),6),4));
                               lines(log(c,10),log(1+c/beta.T.m.baseline,10),lwd=3);
                               lines(log(ell.25,10),lwd=2);
                               lines(log(ell.05,10),lwd=2,lty="dashed");
                               ## points(median.out[1],median.out[2],pch=1)
                               ## points(mean.out[1],mean.out[2],pch=20)
                               points(log(mode.out[1],10),log(mode.out[2],10),pch=20)
},
               key.axes=axis(4,at=pretty(zmat)),
               color=colorRampPalette(brewer.pal(9,"Reds")))

## Plot the R0 as a function of Bt0/(Bt0 + c) and p/mu.i, assuming that 1/3 of midgut infections are established
beta.T.m.c <- 10^seq(-5,-1,0.1)
p.mu.i <- 10^seq(2,6,0.1)
R0 <- expand.grid(beta.T.m.c=beta.T.m.c,p.mu.i=p.mu.i)
R0$R0 <- R0$beta.T.m.c * R0$p.mu.i
zmat <- matrix(R0$R0,nrow=length(beta.T.m.c),ncol=length(p.mu.i))
temp.df <- data.frame(sample.out$beta.m/(sample.out$beta.m+sample.out$c.m),
                      sample.out$p.m/sample.out$mu.m)
median.out <- sapply(temp.df,median)
mean.out <- colMeans(temp.df)
mode.out <- c(MAP.out["beta.m"]/(MAP.out["beta.m"]+MAP.out["c.m"]),
              MAP.out["p.m"]/MAP.out["mu.m"])
cov.out <- cov(temp.df)
ell.25 <- ellipse(mean.out,cov.out,lwd=2,alpha=0.25,draw=FALSE)
ell.05 <- ellipse(mean.out,cov.out,lwd=2,alpha=0.05,draw=FALSE)
## plot
filled.contour(log(beta.T.m.c,10),log(p.mu.i,10),zmat,
               plot.title = title("Midgut R0",
                                  xlab=expression(paste(beta,"/(",beta,"+",c[m],")")),
                                  ylab=expression(paste(p,"/",mu[i]))),
               axes=FALSE,
               plot.axes={axis(1,at=pretty(range(log(beta.T.m.c,10)),6),
                               labels=round(10^pretty(range(log(beta.T.m.c,10)),6),6));
                               axis(2,at=pretty(range(log(p.mu.i,10)),6),
                                    labels=round(10^pretty(range(log(p.mu.i,10)),6),4));
                               lines(log(beta.T.m.c,10),log(1/beta.T.m.c,10),lwd=3);
                               lines(log(ell.25,10),lwd=2);
                               lines(log(ell.05,10),lwd=2,lty="dashed");
                               ## points(median.out[1],median.out[2],col="white",pch=1)
                               ## points(mean.out[1],mean.out[2],col="white",pch=20)
                               points(log(mode.out[1],10),log(mode.out[2],10),pch=20)},
               key.axes=axis(4,at=pretty(zmat)),
               color=colorRampPalette(brewer.pal(9,"Reds")))


## Fit deterministic model
## Run with all four initial conditions and do a weighted sum
## Then run with the IT initial conditions.
## Fit both of these with a poisson likelihood
LL <- function(pars) {
    names(pars) <- names(lower)
    with(as.list(pars),{
        parms <- c(c.b=c.m, beta.b=beta.m,
                   c.m=c.m, beta.m=beta.m, d.m = d.m, p.m=p.m,
                   T0.m=T0.m, mu.m=mu.m, mui.m=mu.m, epsilon.m=epsilon.m,
                   c.s=c.s, beta.s=beta.s, d.s = d.s, p.s=p.s,
                   T0.s=T0.s, mu.s=mu.s, mui.s=mu.s, epsilon.s=epsilon.s)
        
        V.out <- calc.V(state,parms,times)
        V.adj <- (prop.positive.vec - prop.disseminated) * V.out$V.tot + prop.disseminated * V.out$V.tot.pos
        titre.total <- pmax(V.adj[data.hours],0)

        ## Now intrathoracic
        V.tot.it <- calc.V.it(state,parms,times)$V.tot.it
        titre.total.it <- pmax(V.tot.it[data.hours.it],0)
        return(sum(dpois(titre.data,titre.total,log=TRUE),na.rm=TRUE)
               +sum(dpois(titre.data.it,titre.total.it,log=TRUE),na.rm=TRUE))
    })
}

## optim.out <- optim(pars.init,NLL,control=list(maxit=5e3))
## save(optim.out,file="optim_out_barriers.RData")
## load("optim_out_barriers.RData")

## BayesianTools
## lower <- c(c.m=1e-2, beta.m=1e-7, d.m=0.01, p.m=1, T0.m=1, mu.m=1e-2, epsilon.m=0,
##            c.s=1e-2, beta.s=1e-7, d.s=0.01, p.s=1, T0.s=1, mu.s=1e-2, epsilon.s=0)
## upper <- c(c.m=1, beta.m=1e-5, d.m=0.99, p.m=1e4, T0.m=1e4, mu.m=1, epsilon.m=24*2,
##            c.s=1, beta.s=1e-5, d.s=0.99, p.s=1e4, T0.s=1e4, mu.s=1, epsilon.s=24*2)
## bayesianSetup = createBayesianSetup(LL,lower=lower,upper=upper,names=names(lower))
## mcmc.out = runMCMC(bayesianSetup,settings=list(iterations=1e4))

## save(mcmc.out,file="./barrier_mcmc_out.RData")
## load(file="./barrier_mcmc_out.RData")

## Now get new parameters and plot
parms <- MAP(mcmc.out,start=1.5e6)$parametersMAP
parms["c.b"] <- parms["c.m"];parms["beta.b"] <- parms["beta.m"];parms["mui.m"] <- parms["mu.m"];parms["mui.s"] <- parms["mu.s"]
V.out <- calc.V(state,parms,times)
V.adj <- (prop.positive.vec - prop.disseminated) * V.out$V.tot + prop.disseminated * V.out$V.tot.pos
V.tot.it <- calc.V.it(state,parms,times)$V.tot.it
par(mfrow=c(2,1))
plot(times, V.adj,type='l', lwd=3, las=1,xaxs="i",yaxs="i",log="y",bty="n",col="green",
     xlab="time (days)", ylab="titre",ylim=c(min(c(V.adj,titre.data)),max(c(V.out$V.tot.3,titre.data))))
lines(times,V.out$V.tot.1,col="red",lwd=2)
lines(times,V.out$V.tot.2,col="gray",lwd=2)
lines(times,V.out$V.tot.3,col="black",lwd=2)
legend("bottomright",legend=c("Positive midges","Midgut infection barrier","No MIB", "No barriers"),lwd=c(3,2,2,2),
       col=c("green","red","gray","black"),bty="n")
points(times[data.hours],titre.data)
plot(times, V.tot.it,type='l', lwd=2, las=1,xaxs="i",yaxs="i",log="y",bty="n",
     xlab="time (days)", ylab="titre",ylim=c(min(c(V.tot.it,titre.data.it)),max(c(V.tot.it,titre.data.it))))
points(times[data.hours.it],titre.data.it)







## Run stochastic model START HERE
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
plot(-1,0.1,xlim=c(0,14*24),ylim=c(1,1e6),
     xlab="Time (hours)",ylab="Virions",log="y",bty="n",xaxs="i")
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
## trajs.pos[which(trajs[,,c("V.b","V.m","V.l","V.d","V.s")] < 10^0.75)] <- NA
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
plot(-1,0.1,xlim=c(0,14*24),ylim=c(1,1e6),
     xlab="Time (hours)",ylab="Virions",log="y",bty="n",xaxs="i")
trajs.tot <- apply(trajs[,,c("V.b","V.m","V.l","V.d","V.s")],c(1,2),sum)
trajs.totpos <- trajs.tot
trajs.totpos[which(trajs.tot < 10^0.75)] <- NA
for (iii in 1:n.sims) {
    ##print(iii)
    lines(trajs[iii,,"time"],trajs.tot[iii,])
}
lines(trajs[1,,"time"],
      ## rowSums(apply(trajs[,,c("V.b","V.m","V.l","V.d","V.s")],c(2,3),mean)),
      apply(trajs.totpos,2,function(x)mean(x,na.rm=T)),
      col="red",lwd=3)
abline(h=10^0.75,lty="dotted",lwd=2)
abline(h=10^3,lty="dashed",lwd=2)
points(fu$Time.pi..hours,10^fu$log10.mean.titre,col="green",pch=19)
legend("topleft",bty="n",legend=c("Individual midges","Mean of positive midges","Data (Fu, 1999)"),
       lty=c(1,1,NA),col=c("black","red","green"),lwd=c(1,3,NA),pch=c(NA,NA,19))
dev.off()

