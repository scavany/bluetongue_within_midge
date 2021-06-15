## within midge model
wv.BTV.dde = function(t, state, parameters)
{
    with(as.list(c(state, parameters)),{
        if(t < l.s)
            lagval <- 0
        else
            lagval <- lagvalue(t - l.s,1)
        dV.s = -c.s  * V.s + c.s * lagval
        list(c(dV.s))
  })
}

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic$titre <- fu.intrathoracic$titre * log(2) #convert to pfu
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)
V.s0.it <- fu.intrathoracic$titre[1]
V.sinf.it <- mean(tail(fu.intrathoracic$titre,8))
titre.data.it <- floor(fu.intrathoracic$titre+0.5)

state = c(V.s = V.s0.it)
parms = c(c.s = 1e-2)
pars.init <- log(parms)
l.s <- 5

NLL <- function(pars) {
    parms <- exp(pars)
    out <- as.data.frame(dede(state,seq(0,240,1),wv.BTV.dde,parms))
    titre.model.it <- out$V.s
    model.indices <- which(out$time %in% fu.intrathoracic.hours)
    if (any(titre.model.it < 0))
        return(Inf)
    else
        return(-sum(dpois(titre.data.it,titre.model.it[model.indices],log=TRUE)))
}
optim.out <- optimize(NLL,interval=c(-2,6))

parms <- exp(optim.out$par)
out <- as.data.frame(dede(state,seq(0,240,1),wv.BTV.dde,parms))
titre.model.it <- out$V.s

plot(out$time,titre.model.it,type='l',log='y',ylim=range(titre.data.it))
points(fu.intrathoracic.hours,titre.data.it)
