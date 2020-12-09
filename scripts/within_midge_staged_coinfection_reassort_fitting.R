## set working directory
setwd('~/bluetongue_project/aim2/scripts/')

## clear existing workspace
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(deSolve, BayesianTools, data.table, mc2d)

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
## reassortment parameters
growth.ratio <- 100^3
withlike <- 0.75
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

### Fit to el Hussein using multinomial
el.hussein.runtimes <- el.hussein[,.(runtime=max(time)),by=gap]

## range of free parameters
growth.ratio.vec <- 10^c(0,10)
withlike.vec <- c(0,1)
pars <- matrix(c(growth.ratio.vec,withlike.vec),nrow=2,byrow=TRUE,
               dimnames=list(c("growth.ratio","withlike"),
                             c("lower","upper")))

LL <- function(par) {
    growth.ratio <- par[1]
    withlike <- par[2]
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
            times.1 <- seq(0,second.intro,0.1)
            times.2 <- seq(second.intro,runtime,0.1)
            state.1 <- state
            state.1["V.mb"] <- 0
            dat.1 = data.frame(ode(y = state.1, times = times.1,
                                   func = wv.BTV.coinfection.reassort, parms = parameters,
                                   prod.fn.a=prod.wrapper.a,
                                   prod.fn.b=prod.wrapper.b,
                                   prod.fn.r=prod.wrapper.r))
            state.2 <- as.numeric(dat.1[nrow(dat.1),-1])
            names(state.2) <- names(state)
            state.2["V.mb"] <- 1
            dat.2 = data.frame(ode(y = state.2, times = times.2,
                                   func = wv.BTV.coinfection.reassort, parms = parameters,
                                   prod.fn.a=prod.wrapper.a,
                                   prod.fn.b=prod.wrapper.b,
                                   prod.fn.r=prod.wrapper.r))
            dat <- rbind(dat.1,dat.2)
        } else {
            times <- seq(0,runtime,0.1)
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

## sample posterior using MCMC
bayesianSetup = createBayesianSetup(LL,lower=pars[,"lower"],upper=pars[,"upper"])
settings = list(message = T)
out = runMCMC(bayesianSetup, settings = settings)
samples = getSample(out)
save(out,file="mcmc_reassort_out.RData")
exp(colMeans(samples))
plot(out)
