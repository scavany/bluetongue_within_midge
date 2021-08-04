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
       adaptivetau#,psych,bbmle
       )

## load within-vector model of BTV infection
source("./within_midge_barrier_based_fn.R")

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
initial.titre <- 10^fu$log10.mean.titre[1] 

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)
initial.titre.it <- fu.intrathoracic$titre[1]

## Set up timing vectors
times = seq(from = 0, to = 14*24, by = 1)
data.hours.it <- which(times %in% fu.intrathoracic.hours)
model.hours.it <- which(fu.intrathoracic.hours %in% times)
titre.data.it <- floor(fu.intrathoracic$titre[model.hours.it]+0.5)
data.hours <- which(times %in% fu$Time.pi..hours.)
model.hours <- which(fu$Time.pi..hours. %in% times)
titre.data <- floor(10^fu$log10.mean.titre[model.hours] + 0.5)

## Barrier probabilities
p.mib=1/3;p.meb=1/2;p.db=1/5

## Prevalence over time
prop.positive.fn <- function(times,prop.disseminated,k,t0) {
    1 - (1 - prop.disseminated)/(1 + exp(-k*(times - t0)))
}
logistic.optim.fn <- function(par) {
    prop.positive <- prop.positive.fn(times,par[1],par[2],par[3])
    return(sum((prop.positive[data.hours] - fu$Detection.rate[model.hours])^2))
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
prop.disseminated <- min(prop.positive.vec)
t.null <- logistic.optim.out$par[3]

## Calculate midgut clearance based on t.null
V.0 <- floor(initial.titre * log(2) + 0.5)
V.null <- 10^0.75 * log(2)
c.m <- log(V.0/V.null) / t.null

## Oral infection - parms and initial conditions
## The equilibrium is p.s * T.s0 / c.s. This should be set such that
## 80% of the time the equilibrium is below V.null, and
## 20% of the time it is the equilibrium given in the data
## These initial parms are for the 20% - for fitting reasons
V.inf <- mean(log(2) * tail(10^fu$log10.mean.titre))

## Fit deterministic model
## Run with all four initial conditions and do a weighted sum
## Then run with the IT initial conditions.
## Fit both of these with a poisson likelihood
calc.V <- function(state,parms,times) {
    ## Oral infection can result in one of:
    ## 1. constrained to midgut 1 - p.meb
    ## 2. escapes midgut p.meb
    ## 1.
    V.tot.1 <- state["V.m"] * exp(-parms["c.m"] * times)

    ## 2.
    dat.2 = wv.BTV.barrier.det(times,state,parms)
    V.tot.2 <- rowSums(dat.2[,c("V.m","V.s")])
        
    ## Get weighted sum
    V.tot <- V.tot.1*(1-p.meb) + V.tot.2*p.meb
    return(list(V.tot.1=V.tot.1,V.tot.2=V.tot.2,
                V.tot=V.tot))
}
NLL <- function(pars) {
    with(as.list(pars),{
        ## print(exp(lT.s0))
        state <- c(V.m=V.0, T.m=exp(lT.s0)/(1+exp(lT.mult)), E.m=0, I.m=0,
                   V.s=0, T.s=exp(lT.s0), E.s=0, I.s=0)
        beta.m <- c.m * log(1/(1 - p.mib)) / as.numeric(state["T.m"]) / (V.0 - V.null)
        parms <- c(c.m=c.m, beta.m=beta.m, epsilon.m=exp(lepsilon.m),
                   p.m=exp(lp.m), p.s=exp(lp.m)*(1+exp(lp.mult)),
                   c.s=exp(lp.m)*(1+exp(lp.mult))*exp(lT.s0)/V.inf, beta.s=beta.m,
                   epsilon.s=exp(lepsilon.s))
        V.tot <- calc.V(state,parms,times)$V.tot
        titre.total <- V.tot[data.hours]
        return(-sum(dpois(titre.data,titre.total,log=TRUE)))
    })
}
LL <- function(par) {
    names(par) <- names(lower)
    with(as.list(par),{
        ## print(exp(lT.s0))
        state <- c(V.m=V.0, T.m=T.s0/T.mult, E.m=0, I.m=0,
                   V.s=0, T.s=T.s0, E.s=0, I.s=0)
        beta.m <- c.m * log(1/(1 - p.mib)) / as.numeric(state["T.m"]) / (V.0 - V.null)
        parms <- c(c.m=c.m, beta.m=beta.m, epsilon.m=epsilon.m,
                   p.m=p.m, p.s=p.m*p.mult,
                   c.s=p.m*p.mult*T.s0/V.inf, beta.s=beta.s,
                   epsilon.s=epsilon.s)
        V.tot <- calc.V(state,parms,times)$V.tot
        titre.total <- V.tot[data.hours]
        return(sum(dpois(titre.data,titre.total,log=TRUE)))
    })
}

## Fit with optim
pars.init <- c(lT.s0=log(1e5),lT.mult=log(9),
               lepsilon.m=log(2),lepsilon.s=log(3),lp.m=log(1),lp.mult=log(9))
optim.out <- optim(pars.init,NLL,control=list(maxit=5e3))
## save(optim.out,file="optim_out_barriers.RData")
## load("optim_out_barriers.RData")

## Fit with Bayesian Tools
lower <- c(T.s0=1e3,T.mult=1,beta.s=1e-7,
           epsilon.m=1,epsilon.s=1,p.m=1e-3,p.mult=1)
upper <- c(T.s0=1e6,T.mult=1e2,beta.s=1e-3,
           epsilon.m=5,epsilon.s=5,p.m=1e3,p.mult=1e3)
bayesianSetup = createBayesianSetup(LL,lower=lower,upper=upper)
mcmc.out = runMCMC(bayesianSetup,settings=list(iterations=1e5))

## Now get new parameters and plot
samples = getSample(mcmc.out,start=25000)
state <- c(V.m=V.0, T.m=mean(samples[,1])/(1+mean(samples[,2])), E.m=0, I.m=0,
           V.s=0, T.s=mean(samples[,1]), E.s=0, I.s=0)
beta.m <- c.m * log(1/(1 - p.mib)) / as.numeric(state["T.m"]) / (V.0 - V.null)
parms <- c(c.m=c.m, beta.m=beta.m, epsilon.m=mean(samples[,3]),
           p.m=mean(samples[,5]), p.s=mean(samples[,5])*(1+mean(samples[,6])),
           c.s=mean(samples[,5])*(1+mean(samples[,6]))*mean(samples[,1])/V.inf, beta.s=beta.m,
           epsilon.s=mean(samples[,4]))
V.tots <- calc.V(state,parms,times)
plot(times,V.tots$V.tot,log="y",type="l",lwd=3,ylim=c(1,1e4))
lines(times,V.tots$V.tot.1,col="red",lwd=2)
points(fu$Time.pi..hours.,10^fu$log10.mean.titre*log(2))
abline(h=10^0.75,lty="dashed")

## Run stochastic model
n.sims <- 1000
trajs <- array(NA,dim=c(n.sims,dim(wv.BTV.barrier.stoch(times,state,parms))),
               dimnames=list(1:n.sims,
                             rownames(wv.BTV.barrier.stoch(times,state,parms)),
                             colnames(wv.BTV.barrier.stoch(times,state,parms))))

## trajs <- list()
for (iii in 1:n.sims) {
    trajs[iii,,] <- as.matrix(wv.BTV.barrier.stoch(times,state,parms))
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

