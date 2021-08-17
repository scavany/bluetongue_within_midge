library(pacman)
p_load(deSolve, BayesianTools)

## within midge model
wv.BTV.koelle.it = function(t, state, parameters) {
    with(as.list(c(pmax(state,0), parameters)),{
        dV.s = lambda.s * H.s * (1 - dnbinom(0,k,mu=max(P.s,0)/H.s)) - (eta.s + beta.s * H.s) * V.s
        dP.s = beta.s * H.s * V.s - b.s * P.s
        list(c(dV.s, dP.s))
    })
}

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic$titre <- fu.intrathoracic$titre * log(2) #convert to pfu
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)
V.s0.it <- fu.intrathoracic$titre[1]
V.sinf.it <- mean(tail(fu.intrathoracic$titre,8))
titre.data.it <- floor(fu.intrathoracic$titre+0.5)
times.it <- seq(0,max(fu.intrathoracic.hours),1)

state = c(V.s = V.s0.it, P.s = 0)
parms = c(lambda.s = 1e5,
          eta.s = 1e0,
          beta.s = 1e-5,
          b.s = 1e-5,
          H.s = 1e3,
          k = 1e0)
pars.init <- log(parms)
out <- ode(state,times.it,wv.BTV.koelle.it,parms)
## plot(out)

NLL.it <- function(pars) {
    parms <- exp(pars)
    ## if (parms["lambda.s"] < parms["b.s"])
    ##     return(Inf)
    ## else {
    ##     parms["eta.s"] <- parms["beta.s"]*(parms["lambda.s"]/parms["b.s"] - 1)
    ## print(parms)
    out <- as.data.frame(ode(state,times.it,wv.BTV.koelle.it,parms))
    titre.model.it <- out$V.s + out$P.s
    model.indices <- which(out$time %in% fu.intrathoracic.hours)
    if (any(titre.model.it < 0))
        return(Inf)
    else
        return(-sum(dpois(titre.data.it,titre.model.it[model.indices],log=TRUE)))
##}
}
optim.out <- optim(pars.init,NLL.it,control=list(maxit=1e4))

parms.it <- exp(optim.out$par)
out <- as.data.frame(ode(state,times.it,wv.BTV.koelle.it,parms.it))
titre.model.it <- out$V.s + out$P.s

plot(out$time,titre.model.it,type='l',log='y')
lines(out$time,out$V.s,col="red")
lines(out$time,out$P.s,col="green")
points(fu.intrathoracic.hours,titre.data.it)

### Two stage
wv.BTV.koelle = function(t, state, parameters)
{
    with(as.list(c(pmax(state,0), parameters)),{
        dV.m = lambda.m * H.m * (1 - m) * (1 - dnbinom(0,k,mu=max(P.m,0)/H.m)) - (eta.m + beta.m * H.m) * V.m
        dP.m = beta.m * H.m * V.m - b.m * P.m
        dV.s = lambda.m * H.m * m * (1 - dnbinom(0,k,mu=max(P.m,0)/H.m)) +
            lambda.s * H.s * (1 - dnbinom(0,k,mu=max(P.s,0)/H.s)) - (eta.s + beta.s * H.s) * V.s
        dP.s = beta.s * H.s * V.s - b.s * P.s
        list(c(dV.m, dP.m, dV.s, dP.s))
    })
}

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
fu$log10.mean.titre <- log(log(2) * 10^fu$log10.mean.titre,10) # convert to pfu
fu.hours <- fu$Time.pi..hours.
V.m0 <- 10^fu$log10.mean.titre[1]
V.s0 <- 0
V.sinf <- mean(tail(10^fu$log10.mean.titre,10))
titre.data <- floor(10^fu$log10.mean.titre + 0.5)
times <- seq(0,max(fu.hours),1)

state = c(V.m = V.m0, P.m = 0,
          V.s = V.s0, P.s = 0)
state.it = c(V.s = V.s0.it, P.s = 0)
parms = c(lambda.m = as.numeric(parms.it['lambda.s'])/1e6,
          eta.m = as.numeric(parms.it['eta.s']),
          beta.m = as.numeric(parms.it['beta.s']),
          b.m = as.numeric(parms.it['b.s']),
          H.m = as.numeric(parms.it['H.s'])/10,
          parms.it,
          m = 1e-6)
pars.init <- c(log(parms[-length(parms)]), log(1 / (1 / parms[length(parms)] - 1)))

prop.positive.fn <- function(times,prop.disseminated,k,t0) {
    1 - (1 - prop.disseminated)/(1 + exp(-k*(times - t0)))
}
logistic.optim.fn <- function(par) {
    prop.positive <- prop.positive.fn(times,par[1],par[2],par[3])
    return(sum((prop.positive[data.hours] - fu$Detection.rate[model.hours])^2))
}
data.hours <- which(times %in% fu$Time.pi..hours.)
model.hours <- which(fu$Time.pi..hours. %in% times)
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
prop.failure <- prop.positive.vec - prop.disseminated

NLL <- function(pars) {
    parms <- c(exp(pars[-length(pars)]), 1 / (1 + exp(-pars[length(pars)])))
    ## parms["b.m"] <- parms["b.s"]
    ## parms["H.m"] <- parms["H.s"]/10
    out.it <- as.data.frame(ode(state.it,times.it,wv.BTV.koelle.it,parms))
    titre.model.it <- pmax(0,out.it$V.s + out.it$P.s)
    model.indices.it <- which(out.it$time %in% fu.intrathoracic.hours)
    out <- as.data.frame(ode(state,times,wv.BTV.koelle,parms))
    titre.disseminated <- pmax(0,out$V.m + out$P.m  + out$V.s + out$P.s)
    parms["m"] <- 0
    out.fail <- as.data.frame(ode(state,times,wv.BTV.koelle,parms))
    titre.fail <- pmax(out.fail$V.m + out.fail$P.m  + out.fail$V.s + out.fail$P.s,0)
    titre.model <- (titre.disseminated * prop.disseminated + titre.fail * prop.failure) / prop.positive.vec
    model.indices <- which(out$time %in% fu.hours)
    ## if (any(titre.disseminated < 0) | any(titre.model.it < 0) | any(titre.fail < 0))
    ##     return(Inf)
    ## else
    ## test <- sum(dpois(titre.data.it,
    ##                   titre.model.it[model.indices.it],
    ##                   log=TRUE))
    ## if (is.infinite(test) & test > 0) {
    ##     print("inf")
    ##     print(parms)
    ## } else if (is.nan(test)) {
    ##     print("nan")
    ##     print(parms)
    ## }
    return(-sum(dpois(titre.data.it,
                      titre.model.it[model.indices.it],
                      log=TRUE))
           -sum(dpois(titre.data,
                      titre.model[model.indices],
                      log=TRUE)))
}
LL <- function(par) {
    names(par) <- names(parms)
    return(-NLL(par))
}
optim.out <- optim(pars.init,NLL,control=list(maxit=1e4))

## Try MCMC
lower <- c(rep(log(1e-8), 11), log(1 / (1 / 1e-8 - 1)))
upper <- c(rep(log(1e4), 11), log(1 / (1 / (1 - 1e-4) - 1)))
names(lower) <- names(parms)
names(upper) <- names(parms)
bayesianSetup = createBayesianSetup(LL,lower=lower,upper=upper)
mcmc.out = runMCMC(bayesianSetup)
plot(mcmc.out,start=1000)

## Look at results
parms <- c(exp(optim.out$par[-length(optim.out$par)]), 1 / (1 + exp(-optim.out$par[length(optim.out$par)])))
##parms["b.m"] <- parms["b.s"]
##parms["H.m"] <- parms["H.s"]/10
out.it <- as.data.frame(ode(c(V.s=1e5,P.s=0),times,wv.BTV.koelle.it,parms))
titre.model.it <- out.it$V.s + out.it$P.s
out <- as.data.frame(ode(state,times,wv.BTV.koelle,parms))
titre.disseminated <- out$V.m + out$P.m  + out$V.s + out$P.s
parms["m"] <- 0
out.fail <- as.data.frame(ode(state,times,wv.BTV.koelle,parms))
titre.fail <- out.fail$V.m + out.fail$P.m  + out.fail$V.s + out.fail$P.s
titre.model <- (titre.disseminated * prop.disseminated + titre.fail * prop.failure) / prop.positive.vec

par(mfrow=c(2,1))
plot(out.it$time,titre.model.it,type='l',log='y')
points(fu.intrathoracic.hours,titre.data.it)
plot(out$time,titre.model,type='l',log='y')
lines(out$time,titre.disseminated,col="green")
lines(out$time,titre.fail,col="red")
points(fu.hours,titre.data)

### Coinfection model
