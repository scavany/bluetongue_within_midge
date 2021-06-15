## within midge model
wv.BTV = function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    dV.m = -c.m * V.m
    dT.m = -beta.m * T.m * V.m
    dI.m0 = beta.m * T.m * V.m - n * I.m0 / epsilon
    dI.m <- vector(mode="numeric",length=n-1)
    if (n > 1) {
        for (ii in 1:(n-1)) {
            diffval <- get(paste0("I.m",ii-1)) - get(paste0("I.m",ii))
            dI.m[ii] <- n * diffval / epsilon
        }
    }
    dI.mn = n * get(paste0("I.m",n-1)) / epsilon
    dV.h = p.m * I.mn + k * p.s * I.s - c.s * V.h 
    dT.s = -beta.s * T.s * V.h
    dI.s <- beta.s * T.s * V.h
    list(c(dV.m, dT.m, dI.m0, dI.m, dI.mn,
           dV.h, dT.s, dI.s))
  })
}

## within midge model intrathoracic (no eclipse)
wv.BTV.intrathoracic = function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    dV.h = p.s * I.s - c.s * V.h # all infected cells produce virus
    dT.s = -beta.s * T.s * V.h
    dI.s = beta.s * T.s * V.h 
    list(c(dV.h, dT.s, dI.s))
  })
}

wv.BTV.intrathoracic.transitions <- list(
    c(V.h = 1), # Viral production
    c(V.h = -1), # Viral death
    c(T.s = -1, I.s = 1) # Cell infection
)

wv.BTV.intrathoracic.ratefunc <- function(state,parameters,t) {
    with(as.list(c(state, parameters)),{
        return(c(
            p.s * I.s,
            c.s * V.h,
            beta.s * T.s * V.h
        ))
    })
}
