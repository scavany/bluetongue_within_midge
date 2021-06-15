## within midge model, it
wv.BTV.midgut = function(t,c.m,V.m0=1)
{
    return(V.m0*exp(-c.m*t))
}

wv.BTV.secondary = function(t,c.s,V.s0=1,V.sinf=10)
{
    return(pmax(V.s0,V.sinf - (V.sinf - V.s0)*exp(-c.s*t)))
}

wv.BTV.secondary.logistic = function(t,state,parms) {
    with(as.list(c(state,parms)), {
        dV.m <- -c.m*V.m
        dV.s <- g.s * V.s * (1 - V.s / K.s) + alpha * V.m
        list(c(dV.m,dV.s))
  })
}



d1passMIB = function(t,sigma,V.m0,c.m)
{
    return(sigma*V.m0*exp(-c.m*t)*exp(sigma*V.m0/c.m * (exp(-c.m*t) - 1)))
}

d1passMIB.lag = function(t,sigma,V.m0,c.m,lag)
{
    return(ifelse(t<lag,0,sigma*V.m0*exp(-c.m*(t-lag))*exp(sigma*V.m0/c.m * (exp(-c.m*(t-lag)) - 1))))
}

p1passMIB = function(t,sigma,V.m0,c.m)
{
    return(1 - exp(sigma*V.m0/c.m * (exp(-c.m*t) - 1)))
}

p1passMIB.lag = function(t,sigma,V.m0,c.m,lag)
{
    return(ifelse(t<lag,0,1 - exp(sigma*V.m0/c.m * (exp(-c.m*(t-lag)) - 1))))
}

p1passMIB.inverse = function(p,sigma,V.m0,c.m)
{
    return(-log(1+c.m/(V.m0*sigma) * log(1-p))/c.m)
}

r1passMIB = function(n,sigma,V.m0,c.m) {
    x <- runif(n)
    return(p1passMIB.inverse(punif(x),sigma,V.m0,c.m))
}

secondary.tissue.dynamics <- function(t,state,parms) {
    with(as.list(c(state,parms)), {
        prop.A <- (1 - dnbinom(0,size=k,mu=V.sa/C0))/(1-dnbinom(0,size=k,mu=(V.sa + V.sb + V.sr)/C0))
        prop.B <- (1 - dnbinom(0,size=k,mu=V.sb/C0))/(1-dnbinom(0,size=k,mu=(V.sa + V.sb + V.sr)/C0))
        prop.C <- pmax(0,1 - prop.A - prop.B)
        dV.sa <- p.s * prop.A - c.s * V.sa
        dV.sb <- p.s * prop.B - c.s * V.sb
        dV.sr <- p.s * prop.C - c.s * V.sr
        list(c(dV.sa, dV.sb, dV.sr))
  })
}

