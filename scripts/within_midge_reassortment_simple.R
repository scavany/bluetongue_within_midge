library(deSolve)

## Simple virus model
## Vt <- function(t,V0,c,K,r) V0*exp(-c*t) + K/(1+(K-1)*exp(-r*t))
## Vt.inf <- function(t,K,r) K/(1+(K-1)*exp(-r*t))
## t <- seq(0,20,0.1)
## plot(t,Vt(t,1000,1,1e3,1),log="y",type='l')
## V1 <- Vt(t,900,1,1e4,2)
## V2 <- Vt(t,1100,1,1e3,2)
## plot(t,V1,log="y",type='l')
## lines(t,V2,col="red")

## Coinfection parms
simulate_reassortant_ode <- function (V1.vec,V2.vec,t.vec,t02=0, #Vector time-series for V1 and 2 (all same length),
                                      ## times defined at, importation time of 2
                                      N=1e3, k=Inf, p=8, f=1, c=0.7,#N cells, dispersion in MOI, reassortant prod, fitness
                                      dominance=0, #amount first virus dominates second
                                      tstart = 0, tfinal = 20*24, dt = 1)
{
    V2.vec.new <- c(rep(0,sum(t.vec<t02)),V2.vec)
    V1 <- approxfun(t.vec,V1.vec,rule=1:2)
    V2 <- approxfun(c(t.vec,seq(max(t.vec)+mean(diff(t.vec)),max(t.vec)+t02,mean(diff(t.vec)))),
                    V2.vec.new,rule=1:2)
    reassortant_ode <- function(t, y, parms) {
        with(as.list(y, parms),{
            C = (1-dnbinom(0,size=k,mu=V1(t)/N)) * (1-dnbinom(0,size=k,mu=V2(t)/N)) * dnbinom(0,size=k,mu=R/N) + # 1 and 2
                (1-dnbinom(0,size=k,mu=V1(t)/N)) * dnbinom(0,size=k,mu=V2(t)/N) * (1-dnbinom(0,size=k,mu=R/N)) + # 1 and R
                dnbinom(0,size=k,mu=V1(t)/N) * (1-dnbinom(0,size=k,mu=V2(t)/N)) * (1-dnbinom(0,size=k,mu=R/N)) + # 2 and R
                (1-dnbinom(0,size=k,mu=V1(t)/N)) * (1-dnbinom(0,size=k,mu=V2(t)/N)) * (1-dnbinom(0,size=k,mu=R/N)) # all 3
            dRdt = p * C * N - c * R / f
            list(c(dRdt))
        })
    }
    Y0 = c(R = 0)
    timevec = seq(tstart, tfinal, by = dt)
    pars = c(N = N, k = k, p = p, f = f)
    odeout = dede(y = Y0, times = timevec, func = reassortant_ode, 
                  parms = pars, atol = 1e-12, rtol = 1e-12)
    result <- list()
    result$ts <- as.data.frame(odeout)
    result$ts$V1 <- V1(timevec)
    result$ts$V2 <- V2(timevec)
    return(result)
}

out <- simulate_reassortant_ode(V.out$V.tot.pos,V.out$V.tot.pos,times,t02=24)
head(out)
plot(out$ts$time,out$ts$V1, type='l', col="black", lwd=2, bty="n",
     xaxs="i",yaxs="i",log="y",las=1,ylim=c(1,max(out$ts[,2:ncol(out$ts)])))
lines(out$ts$time,out$ts$V2,col="gray",lwd=2,lty=2)
lines(out$ts$time,out$ts$R,col="red",lwd=2,lty=2)
plot(out$ts$time,out$ts$R/(out$ts$V1+out$ts$V2+out$ts$R),type='l',lwd=2)

C = (1-dnbinom(0,size=k,mu=V1(out$ts$time)/N)) * (1-dnbinom(0,size=k,mu=V2(out$ts$time)/N)) * dnbinom(0,size=k,mu=out$ts$R/N) + # 1 and 2
    (1-dnbinom(0,size=k,mu=V1(out$ts$time)/N)) * dnbinom(0,size=k,mu=V2(out$ts$time)/N) * (1-dnbinom(0,size=k,mu=out$ts$R/N)) + # 1 and R
    dnbinom(0,size=k,mu=V1(out$ts$time)/N) * (1-dnbinom(0,size=k,mu=V2(out$ts$time)/N)) * (1-dnbinom(0,size=k,mu=out$ts$R/N)) + # 2 and R
    (1-dnbinom(0,size=k,mu=V1(out$ts$time)/N)) * (1-dnbinom(0,size=k,mu=V2(out$ts$time)/N)) * (1-dnbinom(0,size=k,mu=out$ts$R/N)) # all 3

## try with variable initial conditions and infection times
out <- simulate_reassortant_ode(p=1,t02=2)
head(out)
plot(out$ts$time,out$ts$V1, type='l', col="black", lwd=2, bty="n",
     xaxs="i",yaxs="i",log="y",las=1,ylim=c(1,max(out$ts[,2:ncol(out$ts)])))
lines(out$ts$time,out$ts$V2,col="gray",lwd=2,lty=2)
lines(out$ts$time,out$ts$R,col="red",lwd=2,lty=2)
