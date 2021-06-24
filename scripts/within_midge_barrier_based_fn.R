## hybrid model
wv.BTV.barrier = function(t,
                          N.m=0,N.l=0,N.d=0,N.s=0,N.o=0,
                          p.mib=1,p.meb=1,p.db=1,p.sgib=1,p.sgeb=1,p.totb=1,
                          lambda.m,lambda.l,lambda.d,lambda.s,lambda.o,
                          c.m,c.h,
                          lp.m,lp.l,lp.d,lp.s,lp.o,
                          dp.m,dp.l,dp.d,dp.s,dp.o)
{
   0
}

## deterministic model
wv.BTV.barrier.det = function(t, state, parameters)
{
    with(as.list(c(state, parameters)),{
        dV.b = -c.b * V.b
        dV.m = (1 - d.m) * p.m * I.m - c.m * V.m
        dT.m = -(beta.b * V.b + beta.m * V.m) * T.m 
        dI.m = (beta.b * V.b + beta.m * V.m) * T.m
        dV.l = d.m * p.m * I.m + (1 - d.l) * p.l * I.l - c.l * V.l
        dT.l = -beta.l * T.l * V.l 
        dI.l = beta.l * T.l * V.l 
        dV.d = d.l * p.l * I.l + (1 - d.d) * p.d * I.d - c.d * V.d
        dT.d = -beta.d * T.d * V.d 
        dI.d = beta.d * T.d * V.d 
        dV.s = d.d * p.d * I.d + p.s * I.s - c.s * V.s
        dT.s = -beta.s * T.s * V.s 
        dI.s = beta.s * T.s * V.s 
        list(c(dV.m, dT.m, dI.m, dV.l, dT.l, dI.l,
               dV.d, dT.d, dI.d, dV.s, dT.s, dI.s))
  })
}

## stochastic model
wv.BTV.barrier.stoch <- function(t, state, parameters, barrier.probs){

    ## Calculate available target cells
    state[["T.m"]] <- rbinom(1, state[["N.m"]],
                             1-(1-barrier.probs[["p.mib"]])^(1/state[["N.m"]]))
    state[["T.l"]] <- rbinom(1, state[["N.l"]],
                             1-(1-barrier.probs[["p.meb"]])^(1/state[["N.l"]]))
    state[["T.d"]] <- rbinom(1, state[["N.d"]],
                             1-(1-barrier.probs[["p.db"]])^(1/state[["N.d"]]))
    state[["T.s"]] <- rbinom(1, state[["N.s"]],
                             1-(1-barrier.probs[["p.sgib"]])^(1/state[["N.s"]]))
    
    wv.BTV.transitions <- list(
        c(V.b=-1),
        c(V.m=1),
        c(V.m=-1),
        c(T.m=-1,I.m=1),
        c(V.l=1),
        c(V.l=-1),
        c(T.l=-1,I.l=1),
        c(V.d=1),
        c(V.d=-1),
        c(T.d=-1,I.d=1),
        c(V.s=1),
        c(V.s=-1),
        c(T.s=-1,I.s=1)
    )

    wv.BTV.rates <- function(state, parameters, t)
    {
        with(as.list(c(state, parameters)),{
            return(c(c.b * V.b,
                     (1 - d.m) * p.m * I.m,
                     c.m * V.m,
                     (beta.b * V.b + beta.m * V.m) * T.m,
                     d.m * p.m * I.m + (1 - d.l) * p.l * I.l,
                     c.l * V.l,
                     beta.l * T.l * V.l,
                     d.l * p.l * I.l + (1 - d.d) * p.d * I.d,
                     c.d * V.d,
                     beta.d * T.d * V.d,
                     d.d * p.d * I.d + p.s * I.s,
                     c.s * V.s,
                     beta.s * T.s * V.s
                     ))
        })
    }

    return(simulateModelStochastic(parameters, state, t, wv.BTV.transitions, wv.BTV.rates))
}
