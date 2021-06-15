## set working directory
setwd('~/bluetongue_project/aim2/scripts/')

## clear existing workspace
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(deSolve, BayesianTools, data.table, mc2d,psych)

## load within-vector model of BTV infection
source("./general_functions.R")
source("./btv_culicoides_functions.R")
source("./within_midge_staged_fn.R")
load("./optim_out.RData")
load("./optim_out_intrathoracic.RData")

## Mellor "data"
## mellor <- read.csv("../data/Mellor_points.csv",header=F)
## names(mellor) <- c("days","titre")
## mellor$days <- floor(mellor$days+0.5)
## initial.titre <- 10^mellor$titre[1]

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
initial.titre <- 10^fu$log10.mean.titre[1]
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)

## Samal data
samal <- fread("../data/samal_data.csv")
samal$time <- samal$time*24

## el hussein data
el.hussein <- fread("../data/elHussein_data.csv")
el.hussein$time <- el.hussein$time*24
el.hussein$gap <- el.hussein$gap*24

## fixed parameters
pars <- as.numeric(optim.out$par)
pars.intrathoracic <- as.numeric(optim.out.intrathoracic$par)
n <- 10
beta.s <- exp(pars.intrathoracic[1])/24
p.s <- exp(pars.intrathoracic[2])/24
c.s <- p.s / geometric.mean(tail(fu.intrathoracic$titre))
beta.m <- exp(pars[1])
p.m <- exp(pars[2])
c.m <- exp(pars[3])
epsilon <- exp(pars[4])
tau <- 5
k <- geometric.mean(tail(10^fu$log10.mean.titre))/geometric.mean(tail(fu.intrathoracic$titre)) -
    p.m*(1-exp(-beta.m*initial.titre))/p.s
## reassortment parameters
growth.ratio <- 100^3
withlike <- 0.75

## run model
times = seq(from = 0, to = 20*24, by = 1)
parameters = c(c.m = c.m, c.s = c.s, beta.m = beta.m, beta.s = beta.s, n = n,
               p.m = p.m, p.s = p.s, epsilon = epsilon, tau = tau, k = k)
V.ma <- initial.titre/2; V.mb <- initial.titre/2; T.m <- 1; T.s <- 1 # equations are proportions
state = c(V.ma, V.mb,
          T.m, rep(0,parameters["n"]^2 + 4*parameters["n"] + 1),
          rep(0,3),
          T.s, rep(0,9))
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
                 "I.sbr00")
names(state) <- state.names

### Fit to el Hussein using multinomial
el.hussein.runtimes <- el.hussein[,.(runtime=max(time)),by=gap]

## range of free parameters
growth.ratio.vec <- 10^c(0,10)
withlike.vec <- c(0,1)
tau.vec <- c(1,240)
pars <- matrix(c(growth.ratio.vec,withlike.vec,tau.vec),nrow=3,byrow=TRUE,
               dimnames=list(c("growth.ratio","withlike","tau"),
                             c("lower","upper")))

LL <- function(par) {
    growth.ratio <- par[1]
    withlike <- par[2]
    parameters["tau"] = par[3]
    prod.wrapper.a <- function(i,j,n,r = growth.ratio, wl = withlike){
        production.fn.a(i,j,n,r,wl)
    }
    prod.wrapper.b <- function(i,j,n,r = growth.ratio, wl = withlike){
        production.fn.b(i,j,n,r,wl)
    }
    prod.wrapper.r <- function(i,j,n,r = growth.ratio, wl = withlike){
        production.fn.r(i,j,n,r,wl)
    }
    loglik <- 0
    for (second.intro in unique(el.hussein.runtimes$gap)){
        runtime <- el.hussein.runtimes[gap==second.intro,runtime]
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
            state.2 <- as.numeric(dat.1[nrow(dat.1),-1])
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
        ## Calculate likelihood
        data <- el.hussein[gap==second.intro, .(time,N_17,N_10,N_reassortant,total)]
        dat[,total.virus:=V.ma+V.mb+V.sa+V.sb+V.sr]
        virus.prop <- dat[,.(time,
                             prop.a=(V.ma+V.sa)/total.virus,
                             prop.b=(V.mb+V.sb)/total.virus,
                             prop.r=(V.sr)/total.virus)]
        for (ii in 1: nrow(data)) {
            tt <- as.numeric(data[ii,time])
            loglik <- loglik + dmultinom(as.numeric(data[ii,.(N_17,N_10,N_reassortant)]),
                                         as.numeric(data[ii,.(total)]),
                                         prob=as.numeric(virus.prop[time==tt,
                                                                    .(prop.a,prop.b,prop.r)]),
                                         log=TRUE)
        }
    }
    return(loglik)
}
NLL <- function(par) {-LL(par)}

## sample posterior using MCMC
bayesianSetup = createBayesianSetup(LL,lower=pars[,"lower"],upper=pars[,"upper"])
settings = list(message = T)
out = runMCMC(bayesianSetup, settings = settings)
samples = getSample(out)
save(out,file="mcmc_reassort_out.RData")

### assess MCMC
load("mcmc_reassort_out.RData")
samples = getSample(out,start=1250)
plot(out)
