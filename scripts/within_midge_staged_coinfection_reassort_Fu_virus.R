## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspacer
rm(list = ls())

## set seed
args = c(24,1)
set.seed(args[2])

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(BayesianTools,psych,bbmle,deSolve,grDevices,data.table)

## load within-vector model of BTV infection
source("./within_midge_staged_virus_only_fn.R")
load("virus_only_model_basic.RData")
load("optim_out_virus.RData")
c.s <- as.numeric(exp(optim.out$par['lc.s']))

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
fu$log10.mean.titre <- log(log(2) * 10^fu$log10.mean.titre,10) # convert to pfu
fu.hours <- fu$Time.pi..hours.
V.m0 <- 10^fu$log10.mean.titre[1]
V.s0 <- 0
V.sinf <- mean(tail(10^fu$log10.mean.titre,10))
titre.data <- floor(10^fu$log10.mean.titre + 0.5)

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic$titre <- fu.intrathoracic$titre * log(2) #convert to pfu
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)
V.s0.it <- fu.intrathoracic$titre[1]
V.sinf.it <- mean(tail(fu.intrathoracic$titre,8))
titre.data.it <- floor(fu.intrathoracic$titre+0.5)

### Now for the coinfection
## Load data
## Samal data
samal <- fread("../data/samal_data.csv")
samal$time <- samal$time*24

## el hussein data
el.hussein <- fread("../data/elHussein_data.csv")
el.hussein$time <- el.hussein$time*24
el.hussein$gap <- el.hussein$gap*24

## Pick trajectories for the two viruses
sample.A <- 1:(n/2)
sample.B <- (n/2+1):n
titre.total.A <- titre.total[,sample.A]
titre.total.B <- titre.total[,sample.B]
waits.A <- waits[sample.A]
waits.B <- waits[sample.B]
passes.A <- passes[sample.A]
passes.B <- passes[sample.B]
equilibriums.A <- equilibriums[sample.A]
equilibriums.B <- equilibriums[sample.B]

## work out total combined titre, given lag in midgut
lag = args[1]
titre.midgut.combined <- titre.midgut  + c(rep(0,lag),titre.midgut[1:(length(titre.midgut)-lag)])

## Work out proportion of cells with 0 infections, all of 1 type, and infected with two strains,
## at time of breakthrough
C0 <- 1e3 # number of cells in a single midge
C0.midgut <- C0/10
max.virions <- 1e6
k <- Inf # dispersion in neg bin (set to Inf for poisson)
## prop.0 <- dnbinom(0,size=k,mu=titre.midgut.combined[floor(waits.A+0.5) + 1]/C0.midgut) # uses the fact that timestep is an hour
## prop.A <- rep(0,length(waits.A))
## prop.B <- prop.A
## for (ii in 1:length(prop.A)) {
##     print(ii)
##     prop.A[ii] <- sum(dnbinom(1:max.virions,size=k,mu=titre.midgut.combined[floor(waits.A[ii]+0.5) + 1]/C0.midgut)
##                       * (titre.midgut[floor(waits.A[ii]+0.5) + 1]/titre.midgut.combined[floor(waits.A[ii]+0.5) + 1])^(1:max.virions))
##     prop.B[ii] <- sum(dnbinom(1:max.virions,size=k,mu=titre.midgut.combined[floor(waits.A[ii]+0.5) + 1]/C0.midgut)
##                       * (1 - titre.midgut[floor(waits.A[ii]+0.5) + 1]/titre.midgut.combined[floor(waits.A[ii]+0.5) + 1])^(1:max.virions))
## }
## prop.A <- prop.A/(1-prop.0)
## prop.B <- prop.B/(1-prop.0)
## prop.AB <- 1 - prop.A - prop.B

## estimate the secondary tissue dynamics
times <- seq(1,max(samal$time) - floor(waits.A[args[2]] + 0.5),1)
parms <- c(max.virions=max.virions,k=k,C0=C0,p.s=equilibriums.A[args[2]]*c.s,c.s=c.s)
state <- c(V.sa=titre.midgut[floor(waits.A[args[2]]+0.5) + 1],
           V.sb=(titre.midgut.combined-titre.midgut)[floor(waits.A[args[2]]+0.5) + 1],
           V.sr=0)

## Check if it's all A or B; if so don't run, just set proportion to 0
if (state["V.sb"] == 0 | !passes.A[args[2]]) {
    prop.r <- rep(0,length(times))
    out <- NULL
## Check if A passes. If so, run and get the proportion reassortment
} else  {
    out <- as.data.frame(ode(state,times,secondary.tissue.dynamics,parms))
    prop.r <- out$V.sr / (out$V.sa + out$V.sb + out$V.sr)
} 

## ## Save output
save(prop.r,out, file=paste0("coinfection_output_",args[2],".RData"))

## Fit the k (dispersion) and H (cell population size) params

## plot it all, including average titres
out.list <- list()
prop.r.list <- list()
for (ii in 1:500) {
    load(paste0("coinfection_out/coinfection_output_",ii,".RData"))
    if (!is.null(out)) {
        out.list[[ii]] <- out
    } else {
        out.list[[ii]] <- NA
    }
    prop.r.list[[ii]] <- prop.r
}

tiff("../figures/logistic_positivity.tif",
     res=600,
     compression="lzw",height=600*6,width=600*6)
plot(el.hussein$time[el.hussein$gap==24],el.hussein$prop_reassortant[el.hussein$gap==24],
     ylim=c(0,1),xlim=c(0,384),bty="n",xaxs="i",yaxs="i",las=1,
     xlab="Time (hours)",ylab="Proportion reassortant")
## points(samal$time,samal$prop_reassortant,col="red")
for (ii in 1:500) {
    lines(1:384,c(rep(0,floor(waits.A[ii]+0.5)),
                  prop.r.list[[ii]]))
}
abline(v=24,lty="dashed")
sum(unlist(lapply(prop.r.list,max)) > 0)/length(prop.r.list)
dev.off()
