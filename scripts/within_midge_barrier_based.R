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
       BayesianTools#,psych,bbmle
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
p.mib=1/3;p.meb=1/2;p.db=1/5
p.sgib=1;p.sgeb=1;p.totb=0
barrier.probs <- c(p.mib=p.mib,p.meb=p.meb,p.db=p.db,
                   p.sgib=p.sgib,p.sgeb=p.sgeb,p.totb=p.totb)

## Oral infection - parms and initial conditions
c.m <- 0.2; c.l <- 20; c.s <- 0.05; epsilon <- 2*24
beta <- 1e-5; mu = 1/(3*24); T0 <- 1000
d <- 0.5; p.m <- 10; p.s <- 100
parms <- c(c.b=c.m, beta.b=beta,
           c.m=c.m, epsilon.m=epsilon, d.m = d, p.m=p.m, beta.m=beta, T0.m=T0, mu.m=mu, mui.m=mu,
           c.l=c.l, epsilon.l=epsilon, d.l = d, p.l=p.m, beta.l=beta, T0.l=T0, mu.l=mu, mui.l=mu,
           c.d=c.s, epsilon.d=epsilon, d.d = d, p.d=p.s, beta.d=beta, T0.d=T0, mu.d=mu, mui.d=mu,
           c.s=c.s, epsilon.s=epsilon, d.s = d, p.s=p.s, beta.s=beta, T0.s=T0, mu.s=mu, mui.s=mu)
pars.init <- c(lc.m=log(c.m),lc.l=log(c.l),lc.s=log(c.s),lbeta=log(beta),logitd=qlogis(d),
               lp.m=log(p.m),lp.s=log(p.s),lepsilon=log(epsilon),
               lT.0=log(1e3),lmu=log(mu))

state <- c(V.b = floor(initial.titre*log(2) + 0.5),
           V.m=0, N.m=T0, T.m=0, E.m=0, I.m=0,
           V.l=0, N.l=T0, T.l=0, E.l=0, I.l=0,
           V.d=0, N.d=T0, T.d=0, E.d=0, I.d=0, 
           V.s=0, N.s=T0, T.s=0, E.s=0, I.s=0)


## Fit deterministic model
## Run with all four initial conditions and do a weighted sum
## Then run with the IT initial conditions.
## Fit both of these with a poisson likelihood
calc.V <- function(state,parms,times) {
    ## Oral infection can result in one of:
    ## 1. no infection (1 - p.mib)
    ## 2. constrained to midgut p.mib(1 - p.meb)
    ## 3. escapes midgut, but not disseminated p.mib*p.meb*(1-p.db)
    ## 4. disseminated infection p.mib*p.meb*p.db
    ## 1. no infection (1 - p.mib)
    state.1 <- state
    dat.1 = wv.BTV.barrier.det(times,state.1,parms)
    V.tot.1 <- rowSums(dat.1[,c("V.b","V.m","V.l","V.d","V.s")])
    ## 2. constrained to midgut p.mib(1 - p.meb)
    state.2 <- state.1
    state.2["T.m"] <- state.1["N.m"]*(1-(1-p.mib)^(1/state.1["N.m"]))/p.mib
    dat.2 = wv.BTV.barrier.det(times,state.2,parms)
    V.tot.2 <- rowSums(dat.2[,c("V.b","V.m","V.l","V.d","V.s")])
    
    ## 3. escapes midgut, but not disseminated p.mib*p.meb
    ## May need to tune the parameters to get the dissemination barrier, 
    ## as it's a consequence of number of fat body cells productive and local viral clearance
    state.3 <- state.2
    state.3["T.l"] <- state.2["N.l"]#*(1-(1-p.meb)^(1/state.2["N.l"]))/p.meb
    parms.3 <- parms
    parms.3["d.l"] <- 0
    dat.3 = wv.BTV.barrier.det(times,state.3,parms.3)
    V.tot.3 <- rowSums(dat.3[,c("V.b","V.m","V.l","V.d","V.s")])
    
    ## 4. disseminated infection p.mib*p.meb*p.db (written assuming p.sgib=1)
    state.4 <- state.3
    state.4["T.d"] <- state.3["N.d"]#*(1-(1-p.db)^(1/state.3["N.d"]))/p.db
    state.4["T.s"] <- state.3["N.s"]
    dat.4 = wv.BTV.barrier.det(times,state.4,parms)
    V.tot.4 <- rowSums(dat.4[,c("V.b","V.m","V.l","V.d","V.s")])
    
    ## Get weighted sum
    V.tot <- V.tot.1*(1-p.mib) + V.tot.2*p.mib*(1-p.meb) +
        V.tot.1*p.mib*p.meb*(1-p.db) + V.tot.4*p.mib*p.meb*p.db
    return(list(V.tot.1=V.tot.1,V.tot.2=V.tot.2,V.tot.3=V.tot.3,
                V.tot.4=V.tot.4,V.tot=V.tot))
}
calc.V.it <- function(state,parms,times) {
    ## Oral infection can result in one of:
    state.it <- state
    state.it["V.b"] <- 0
    state.it["V.l"] <- initial.titre.it*log(2)/2
    state.it["V.d"] <- initial.titre.it*log(2)/2
    state.it["T.l"] <- state.it["N.l"]*(1-(1-p.meb)^(1/state.it["N.l"]))/p.meb
    state.it["T.d"] <- state.it["N.d"]#*(1-(1-p.db)^(1/state.it["N.d"]))/p.db
    state.it["T.s"] <- state.it["N.s"]

    dat.it = wv.BTV.barrier.det(times,state.it,parms)
    V.tot.it <- rowSums(dat.it[,c("V.b","V.m","V.l","V.d","V.s")])
    return(list(V.tot.it=V.tot.it))
}
NLL <- function(pars) {
    ## print(pars)
    ## Secondary tissues, intrathoracic first
    with(as.list(pars),{
        ## Oral infection first
        c.m <- as.numeric(exp(lc.m))
        c.l <- as.numeric(exp(lc.l))
        c.s <- as.numeric(exp(lc.s))
        beta <- as.numeric(exp(lbeta))
        d <- as.numeric(plogis(logitd))
        p.m <- as.numeric(exp(lp.m))
        p.s <- as.numeric(exp(lp.s))
        T.0 <- as.numeric(exp(lT.0))
        epsilon <- as.numeric(exp(lepsilon))
        mu <- as.numeric(exp(mu))
       
        ## state["N.m"] <- floor(as.numeric(exp(lN))+0.5)
        ## state["N.l"] <- floor(as.numeric(exp(lN))+0.5)
        ## state["N.d"] <- floor(as.numeric(exp(lN))+0.5)
        ## state["N.s"] <- floor(as.numeric(exp(lN.s))+0.5)
        parms <- c(c.b=c.m, beta.b=beta,
                   c.m=c.m, epsilon.m=epsilon, d.m = d, p.m=p.m, beta.m=beta, T0.m=T0, mu.m=mu, mui.m=mu,
                   c.l=c.l, epsilon.l=epsilon, d.l = d, p.l=p.m, beta.l=beta, T0.l=T0, mu.l=mu, mui.l=mu,
                   c.d=c.s, epsilon.d=epsilon, d.d = d, p.d=p.s, beta.d=beta, T0.d=T0, mu.d=mu, mui.d=mu,
                   c.s=c.s, epsilon.s=epsilon, d.s = d, p.s=p.s, beta.s=beta, T0.s=T0, mu.s=mu, mui.s=mu)
        
        V.out <- calc.V(state,parms,times)
        pos.adj <- (prop.positive.vec - prop.disseminated)/(1 - prop.disseminated)
        titre.total <- (V.out$V.tot + V.out$V.tot.1 * p.mib * pos.adj)[data.hours]

        ## Now intrathoracic
        ##V.tot.it <- calc.V.it(state,parms,times)$V.tot.it
        ##titre.total.it <- V.tot.it[data.hours.it]
        return(-sum(dpois(titre.data,titre.total,log=TRUE)))
        ##-sum(dpois(titre.data.it,titre.total.it,log=TRUE)))
    })
}

optim.out <- optim(pars.init,NLL,control=list(maxit=5e3))
save(optim.out,file="optim_out_barriers.RData")
load("optim_out_barriers.RData")

## Now get new parameters and plot
c.m <- as.numeric(exp(optim.out$par["lc.m"]))
c.l <- as.numeric(exp(optim.out$par["lc.l"]))
c.s <- as.numeric(exp(optim.out$par["lc.s"]))
beta <- as.numeric(exp(optim.out$par["lbeta"]))
d <- as.numeric(plogis(optim.out$par["logitd"]))
## d <- as.numeric(exp(optim.out$par["ld"]))
p.m <- as.numeric(exp(optim.out$par["lp.m"]))
p.s <- as.numeric(exp(optim.out$par["lp.s"]))
epsilon <- as.numeric(exp(optim.out$par["lepsilon"]))
T.0 <- as.numeric(exp(optim.out$par["lT.0"]))
mu <- as.numeric(exp(optim.out$par["lmu"]))
parms <- c(c.b=c.m, beta.b=beta,
           c.m=c.m, epsilon.m=epsilon, d.m = d, p.m=p.m, beta.m=beta, T0.m=T0, mu.m=mu, mui.m=mu,
           c.l=c.l, epsilon.l=epsilon, d.l = d, p.l=p.m, beta.l=beta, T0.l=T0, mu.l=mu, mui.l=mu,
           c.d=c.s, epsilon.d=epsilon, d.d = d, p.d=p.s, beta.d=beta, T0.d=T0, mu.d=mu, mui.d=mu,
           c.s=c.s, epsilon.s=epsilon, d.s = d, p.s=p.s, beta.s=beta, T0.s=T0, mu.s=mu, mui.s=mu)
V.tots <- calc.V(state,parms,times)
pos.adj <- (prop.positive.vec - prop.disseminated)/(1 - prop.disseminated)
titre.total <- (V.tots$V.tot + V.tots$V.tot.1 * p.mib * pos.adj)
plot(times,titre.total,log="y",type="l",lwd=3,ylim=c(1,1e10))
lines(times,V.tots$V.tot.1,col="red",lwd=2)
points(fu$Time.pi..hours.,10^fu$log10.mean.titre*log(2))
abline(h=10^0.75,lty="dashed")
##V.tots.it <- calc.V.it(state,parms,times)
##plot(times,V.tots.it$V.tot.it,log="y",type="l",lwd=3,ylim=c(1e2,1e6))
##lines(times,V.tots$V.tot.1,col="red",lwd=2)
##points(fu.intrathoracic.hours,fu.intrathoracic$titre*log(2))
##abline(h=10^0.75,lty="dashed")

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

