library(deSolve)

init.state <- c(G = 1, R = 0, C = 0, E = 2000)

parms <- c(alpha = 5e-1, kappa = 5e-4, gamma = 1e-6, rho = 3e-4,
           delta.g = 1e-2, delta.r = 1e-2, delta.c = 1e-2)

cooperation.odes = function(t, state, parameters)
{
    with(as.list(c(state, parameters)),{
        dG <- 2 * alpha * C - gamma * G * R * E - rho * G * R - delta.g * G
        dR <- alpha * C + kappa * G * E - gamma * G * R * E - rho * G * R - delta.r * R
        dC <- gamma * G * R * E - alpha * C - delta.c * C
        dE <- -E * (gamma * G * R + kappa * G)
        dV <- rho * G * R
        list(c(dG, dR, dC, dE))
  })
}

out <- as.data.frame(ode(init.state,seq(1,2000,0.1),cooperation.odes,parms))

plot(out$time,out$G,col="red",log="y",ylim=c(1,2000),ylab="Abundance",type="l",lwd=3,bty="n",
     xlab="Time")
lines(out$time,out$R,col="green",lwd=3)
lines(out$time,out$C,col="blue",lwd=3)
lines(out$time,out$E,col="orange",lwd=3)

plot(out$time,out$G*out$R*parms["rho"],col="red",ylab="Rate of virion production",
     xlab="Time",type="l",lwd=3,bty="n")

intercellular.odes = function(t, state, parameters)
{
    with(as.list(parameters),{
        S <- as.numeric(state["S"])
        I <- as.numeric(state[paste0("I.",1:Kappa)])
        V <- as.numeric(state["V"])
        dI <- rep(0,Kappa)
        dS <- lambda - beta * S * V - delta.s * S
        dI[1] <- beta * S * V - (delta.i + 1) * I[1]
        dI[2:(Kappa-1)] <- I[1:(Kappa-2)] - (delta.i + 1) * I[2:(Kappa-1)]
        dI[Kappa] <- I[Kappa-1] - delta.i * I[Kappa]
        dV <- sum(Vdot * I) - beta * S * V - delta.v * V
        list(c(dS,dI,dV))
    })
}

parms.inter <- c(Kappa=2e3,beta=3e-6,lambda=25,delta.s=2e-3,delta.i=7e-3,
                 delta.v=4e-2)
Vdot <- parms["rho"] * out$G[out$time %in% 1:parms.inter["Kappa"]] * out$R[out$time %in% 1:parms.inter["Kappa"]]
init.state.inter <- c(10^4,1,rep(0,parms.inter["Kappa"]-1),0)
names(init.state.inter) <- c("S",paste0("I.",1:parms.inter["Kappa"]),"V")

out.inter <- as.data.frame(ode(init.state.inter,seq(1,1e4,0.1),intercellular.odes,parms.inter))

plot(out.inter$time,out.inter$V)
