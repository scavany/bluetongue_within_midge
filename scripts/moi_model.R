library(deSolve)

koelle <- function(t, state, parameters)
{
    with(as.list(c(state, parameters)),{
        dH = b*(H0 - H)
        dV = lambda*P - eta*V - beta*H*V
        dP = beta*H*V - b*P
        list(c(dH, dV, dP))
    })
}

## fixed parameters
k <- 0.89
lambda <- 26.75
eta <- 2.91e3
beta <- 7.5e-5
H0=3.5e11
b <- lambda*beta*H0/(eta+beta*H0)

## run model
times = seq(from = 0, to = 20, by = 0.01)
parameters = c(k=k,lambda=lambda,beta=beta,eta=eta,H0=H0,b=b)
state = c(H=H0,V=10^-2.06,P=0)

out <- ode(y = state, times = times, func = koelle, parms = parameters)

plot(out)



## My model
