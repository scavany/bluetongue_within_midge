## deterministic models
wv.BTV.barrier.det = function(t, state, parameters)
{
    wv.BTV.odes <- function(t,state,parameters) {
        with(as.list(c(state, parameters)),{
            dV.b = -beta.b * V.b * (T.m + E.m) - c.b * V.b

            dV.m = -beta.m * V.m * (T.m + E.m) + d.s * p.s * (I.s + I.ss) +
                (1 - d.m) * p.m * (I.m + I.mm) - c.m * V.m
            dT.m = -(beta.b * V.b + beta.m * V.m) * T.m + mu.m * (T0.m - T.m)
            if (epsilon.m == 0) {
                dE.m = -mu.m * E.m
                dI.m = (beta.b * V.b + beta.m * V.m) * T.m - mui.m * I.m
            } else {
                dE.m = (beta.b * V.b + beta.m * V.m) * (T.m - E.m) - E.m / epsilon.m - mu.m * E.m
                dI.m = E.m / epsilon.m - mui.m * I.m
                dE.mm = (beta.b * V.b + beta.m * V.m) * E.m - E.mm / epsilon.m - mu.m * E.mm
                dI.mm = E.mm / epsilon.m - mui.m * I.mm
            }

            dV.s = -beta.s * (T.s + E.s) * V.s + (1 - d.s) * p.s * (I.s + I.ss) +
                d.m * p.m * (I.m + I.mm) - c.s * V.s
            dT.s = -beta.s * T.s * V.s + mu.s * (T0.s - T.s)
            if (epsilon.s == 0) {
                dE.s = -mu.s * E.s
                dI.s = beta.s * V.s * T.s - mui.s * I.s
            } else {
                dE.s = beta.s * V.s * (T.s - E.s) - E.s / epsilon.s - mu.s * E.s
                dI.s = E.s / epsilon.s - mui.s * I.s
                dE.ss = beta.s * V.s * E.s - E.ss / epsilon.s - mu.s * E.ss
                dI.ss = E.ss / epsilon.s - mui.s * I.ss
            }
            list(c(dV.b,
                   dV.m, dT.m, dE.m, dE.mm, dI.m, dI.mm,
                   dV.s, dT.s, dE.s, dE.ss, dI.s, dI.ss))
        })
    }
    out <- as.data.frame(ode(state,t,wv.BTV.odes,parameters))
    return(out)
}

calc.V <- function(state,parms,times) {
    ## Oral infection can result in one of:
    ## 1. no infection (1 - p.mib)
    ## 2. constrained to midgut p.mib(1 - p.meb*p.db)
    ## 3. Fully disseminated p.mib*p.meb*p.db
    ## 1. no infection (1 - p.mib)
    state.1 <- state
    parms.1 <- parms
    parms.1["T0.m"] <- 0
    parms.1["T0.s"] <- 0
    dat.1 = wv.BTV.barrier.det(times,state.1,parms.1)
    V.tot.1 <- rowSums(dat.1[,grepl("V.",names(dat.1))])

    ## 2. constrained to midgut p.mib(1 - p.meb*p.db)
    state.2 <- state.1
    state.2["T.m"] <- parms[["T0.m"]]
    parms.2 <- parms
    ##parms.2["d.m"] <- 0 ## Could also do this with T0
    parms.2["T0.s"] <- 0
    dat.2 = wv.BTV.barrier.det(times,state.2,parms.2)
    V.tot.2 <- rowSums(dat.2[,grepl("V.",names(dat.2))])
    
    ## 3. disseminated infection p.mib*p.meb*p.db (written assuming p.sgib=1)
    ## May need to tune the parameters to get the dissemination barrier, 
    ## as it's a consequence of number of fat body cells productive and local viral clearance
    state.3 <- state.2
    state.3["T.s"] <- parms[["T0.s"]]
    parms.3 <- parms
    dat.3 = wv.BTV.barrier.det(times,state.3,parms.3)
    V.tot.3 <- rowSums(dat.3[,grepl("V.",names(dat.3))])
        
    ## Get weighted sum
    V.tot <- V.tot.1*(1-p.mib) + V.tot.2*p.mib*(1-p.db) +
        V.tot.3*p.mib*p.db
    V.tot.pos <- V.tot.2*(1-p.db) + V.tot.3*p.db
    return(list(V.tot.1=V.tot.1,V.tot.2=V.tot.2,V.tot.3=V.tot.3,
                V.tot=V.tot,V.tot.pos=V.tot.pos))
}
calc.V.it <- function(state,parms,times) {
    ## Oral infection can result in one of:
    state.it <- state
    state.it["V.b"] <- 0
    state.it["V.s"] <- initial.titre.it
    state.it["T.s"] <- parms[["T0.s"]]

    dat.it = wv.BTV.barrier.det(times,state.it,parms)
    V.tot.it <- rowSums(dat.it[,grepl("V.",names(dat.it))])
    return(list(V.tot.it=V.tot.it))
}

## stochastic models
wv.BTV.barrier.stoch <- function(t, state, parameters, barrier.probs, pass.mib=NA, pass.db=NA){

    ## Calculate available target cells - would make more sense to do this outside of function or in wrapper function?
    if(is.na(pass.mib)) pass.mib <- (runif(1) < barrier.probs[["p.mib"]])
    if(is.na(pass.db)) pass.db <- (runif(1) < barrier.probs[["p.db"]])
    parameters[["T0.m"]] <- ifelse(pass.mib,parameters[["T0.m"]],0)
    state[["T.m"]] <- floor(parameters[["T0.m"]] + 0.5)
    parameters[["T0.s"]] <- ifelse(pass.db,parameters[["T0.s"]],0)
    state[["T.s"]] <- floor(parameters[["T0.s"]] + 0.5)
    ##parameters[["d.m"]] <- ifelse(pass.db,parameters[["d.m"]],0)
    
    ## Transitions matrix
    epsilon.m <- parameters[["epsilon.m"]]
    epsilon.s <- parameters[["epsilon.s"]]
    wv.BTV.transitions <- list(
        c(V.b=-1), 
        c(T.m=-1,V.b=-1,E.m=ifelse(epsilon.m==0,0,1),I.m=ifelse(epsilon.m==0,1,0)),

        c(V.m=1), 
        c(V.m=-1), 
        c(T.m=-1,V.m=-1,E.m=ifelse(epsilon.m==0,0,1),I.m=ifelse(epsilon.m==0,1,0)),
        c(E.m=ifelse(epsilon.m==0,0,-1),I.m=ifelse(epsilon.m==0,0,1)),
        c(T.m=1),
        c(T.m=-1),
        c(E.m=-1),
        c(I.m=-1),
        
        c(V.s=1), 
        c(V.s=-1), 
        c(T.s=-1,V.s=-1,E.s=ifelse(epsilon.s==0,0,1),I.s=ifelse(epsilon.s==0,1,0)),
        c(E.s=ifelse(epsilon.s==0,0,-1),I.s=ifelse(epsilon.s==0,0,1)),
        c(T.s=1),
        c(T.s=-1),
        c(E.s=-1),
        c(I.s=-1)
    )
    ## Rates function
    wv.BTV.rates <- function(state, parameters, t)
    {
        with(as.list(c(state, parameters)),{
            return(c(c.b * V.b,
                     beta.b * V.b * T.m,
                     
                     (1 - d.m) * p.m * I.m + d.s * p.s * I.s,
                     c.m * V.m,
                     beta.m * V.m * T.m,
                     ifelse(epsilon.m==0,0,E.m / epsilon.m),
                     T0.m * mu.m,
                     T.m * mu.m,
                     E.m * mu.m,
                     I.m * mui.m,

                     d.m * p.m * I.m + p.s * (1 - d.s) * I.s,
                     c.s * V.s,
                     beta.s * T.s * V.s,
                     ifelse(epsilon.s==0,0,E.s / epsilon.s),
                     T0.s * mu.s,
                     T.s * mu.s,
                     E.s * mu.s,
                     I.s * mui.s
                     ))
        })
    }
    
    out.temp <- as.data.frame(ssa.adaptivetau(state, wv.BTV.transitions, wv.BTV.rates,
                                              parameters, tf=diff(range(t))))
    out.temp$time <- out.temp$time + min(t)
    out <- cbind(time = t, apply(out.temp[, -1], 2, function(col) {
        approx(x = out.temp[, 1], y = col, xout = t, method = "constant")$y
    }))
    return(as.data.frame(out))
}

wv.BTV.coinfection.reassort <- function(t, state, parameters, withlike=0) {
    wv.BTV.coinfection.reassort.odes <- function(t, state, parameters, withlike=0)
    {
        with(as.list(c(state, parameters)),{
            ## Bloodmeal virus
            dV.ba = -beta.b * V.ba * (T.m + E.ma + E.mb + E.mr) - c.b * V.ba
            dV.bb = -beta.b * V.bb * (T.m + E.ma + E.mb + E.mr) - c.b * V.bb
            dV.br = -beta.b * V.br * (T.m + E.ma + E.mb + E.mr) - c.b * V.br
            
            ## Midgut virus
            dV.ma = (-beta.m * V.ma * (T.m + E.ma + E.mb + E.mr) + d.s * p.s * (I.sa + I.saa) +
                     (1 - d.m) * p.m * (I.ma + I.maa) - c.m * V.ma +
                     rho.1(withlike) * (d.s * p.s * (I.sar + I.sab) +
                                        (1 - d.m) * p.m * (I.mar + I.mab)))
            dV.mb = (-beta.m * V.mb * (T.m + E.ma + E.mb + E.mr) + d.s * p.s * (I.sb + I.sbb) +
                     (1 - d.m) * p.m * (I.mb + I.mbb) - c.m * V.mb +
                     rho.1(withlike) * (d.s * p.s * (I.sbr + I.sab) +
                                        (1 - d.m) * p.m * (I.mbr + I.mab)))
            dV.mr = (-beta.m * V.mr * (T.m + E.ma + E.mb + E.mr) + d.s * p.s * (I.sr + I.srr) + 
                     (1 - d.m) * p.m * (I.mr + I.mrr) - c.m * V.mr +
                     rho.1(withlike) * (d.s * p.s * (I.sar + I.sbr) +
                                        (1 - d.m) * p.m * (I.mar + I.mbr)) +
                     rho.r(withlike) * (d.s * p.s * (I.sar + I.sbr + I.sab) +
                                        (1 - d.m) * p.m * (I.mar + I.mbr + I.mab)))
            
            ## Midgut cells
            dT.m = -(beta.b * (V.ba + V.bb + V.br) +
                     beta.m * (V.ma + V.mb + V.mr)) * T.m + mu.m * (T0.m - T.m)
            
            dE.ma = (beta.b * V.ba + beta.m * V.ma) * T.m - (beta.b * V.bb + beta.m * V.mb) * E.ma -
                (beta.b * V.br + beta.m * V.mr) * E.ma - (beta.b * V.ba + beta.m * V.ma) * E.ma -
                E.ma / epsilon.m - mu.m * E.ma
            dE.mb = (beta.b * V.bb + beta.m * V.mb) * T.m - (beta.b * V.ba + beta.m * V.ma) * E.mb -
                (beta.b * V.br + beta.m * V.mr) * E.mb - (beta.b * V.bb + beta.m * V.mb) * E.mb -
                E.mb / epsilon.m - mu.m * E.mb
            dE.mr = (beta.b * V.br + beta.m * V.mr) * T.m - (beta.b * V.ba + beta.m * V.ma) * E.mr -
                (beta.b * V.bb + beta.m * V.mb) * E.mr - (beta.b * V.br + beta.m * V.mr) * E.mr -
                E.mr / epsilon.m - mu.m * E.mr
            dE.maa = (beta.b * V.ba + beta.m * V.ma) * E.ma - E.maa / epsilon.m - mu.m * E.maa
            dE.mbb = (beta.b * V.bb + beta.m * V.mb) * E.mb - E.mbb / epsilon.m - mu.m * E.mbb
            dE.mrr = (beta.b * V.br + beta.m * V.mr) * E.mr - E.mrr / epsilon.m - mu.m * E.mrr
            dE.mab = (beta.b * V.ba + beta.m * V.ma) * E.mb +
                (beta.b * V.bb + beta.m * V.mb) * E.ma - E.mab / epsilon.m - mu.m * E.mab
            dE.mar = (beta.b * V.ba + beta.m * V.ma) * E.mr +
                (beta.b * V.br + beta.m * V.mr) * E.ma - E.mar / epsilon.m - mu.m * E.mar
            dE.mbr = (beta.b * V.bb + beta.m * V.mb) * E.mr +
                (beta.b * V.br + beta.m * V.mr) * E.mb - E.mbr / epsilon.m - mu.m * E.mbr

            dI.ma = E.ma / epsilon.m - mui.m * I.ma
            dI.mb = E.mb / epsilon.m - mui.m * I.mb
            dI.mr = E.mr / epsilon.m - mui.m * I.mr
            dI.maa = E.maa / epsilon.m - mui.m * I.maa
            dI.mbb = E.mbb / epsilon.m - mui.m * I.mbb
            dI.mrr = E.mrr / epsilon.m - mui.m * I.mrr
            dI.mab = E.mab / epsilon.m - mui.m * I.mab
            dI.mar = E.mar / epsilon.m - mui.m * I.mar
            dI.mbr = E.mbr / epsilon.m - mui.m * I.mbr
            
            ## Secondary tissue virus
            dV.sa = (-beta.s * V.sa * (T.s + E.sa + E.sb + E.sr) + d.m * p.m * (I.ma + I.maa) +
                (1 - d.s) * p.s * (I.sa + I.saa) - c.s * V.sa +
                rho.1(withlike) * (d.m * p.m * (I.mar + I.mab) + (1 - d.s) * p.s * (I.sar + I.sab)))
            dV.sb = (-beta.s * V.sb * (T.s + E.sa + E.sb + E.sr) + d.m * p.m * (I.mb + I.mbb) +
                (1 - d.s) * p.s * (I.sb + I.sbb) - c.s * V.sb +
                rho.1(withlike) * (d.m * p.m * (I.mbr + I.mab) + (1 - d.s) * p.s * (I.sbr + I.sab)))
            dV.sr = (-beta.s * V.sr * (T.s + E.sa + E.sb + E.sr) + d.m * p.m * (I.mr + I.mrr) +
                (1 - d.s) * p.s * (I.sr + I.srr) - c.s * V.sr +
                rho.1(withlike) * (d.m * p.m * (I.mar + I.mbr) + (1 - d.s) * p.s * (I.sar + I.sbr)) +
                rho.r(withlike) * (d.m * p.m * (I.mar + I.mbr + I.mab) +
                                   (1 - d.s) * p.s * (I.sar + I.sbr + I.sab)))
            
            ## Secondary tissue cells
            dT.s = -beta.s * (V.sa + V.sb + V.sr) * T.s + mu.s * (T0.s - T.s)

            dE.sa = beta.s * V.sa * T.s - beta.s * V.sb * E.sa - beta.s * V.sr * E.sa -
                beta.s * V.sa * E.sa - E.sa / epsilon.s - mu.s * E.sa
            dE.sb = beta.s * V.sb * T.s - beta.s * V.sa * E.sb - beta.s * V.sr * E.sb -
                beta.s * V.sb * E.sb - E.sb / epsilon.s - mu.s * E.sb
            dE.sr = beta.s * V.sr * T.s - beta.s * V.sa * E.sr - beta.s * V.sb * E.sr -
                beta.s * V.sr * E.sr - E.sr / epsilon.s - mu.s * E.sr
            dE.saa = beta.s * V.sa * E.sa - E.saa / epsilon.s - mu.s * E.saa
            dE.sbb = beta.s * V.sb * E.sb - E.sbb / epsilon.s - mu.s * E.sbb
            dE.srr = beta.s * V.sr * E.sr - E.srr / epsilon.s - mu.s * E.srr
            dE.sab = beta.s * V.sa * E.sb + beta.s * V.sb * E.sa - E.sab / epsilon.s - mu.s * E.sab
            dE.sar = beta.s * V.sa * E.sr + beta.s * V.sr * E.sa - E.sar / epsilon.s - mu.s * E.sar
            dE.sbr = beta.s * V.sb * E.sr + beta.s * V.sr * E.sb - E.sbr / epsilon.s - mu.s * E.sbr
            
            dI.sa = E.sa / epsilon.s - mui.s * I.sa
            dI.sb = E.sb / epsilon.s - mui.s * I.sb
            dI.sr = E.sr / epsilon.s - mui.s * I.sr
            dI.saa = E.saa / epsilon.s - mui.s * I.saa
            dI.sbb = E.sbb / epsilon.s - mui.s * I.sbb
            dI.srr = E.srr / epsilon.s - mui.s * I.srr
            dI.sab = E.sab / epsilon.s - mui.s * I.sab
            dI.sar = E.sar / epsilon.s - mui.s * I.sar
            dI.sbr = E.sbr / epsilon.s - mui.s * I.sbr
            
            list(c(dV.ba, dV.bb, dV.br,
                   dV.ma, dV.mb, dV.mr,
                   dT.m,
                   dE.ma,dE.mb,dE.mr,dE.maa,dE.mbb,dE.mrr,dE.mab,dE.mar,dE.mbr,
                   dI.ma,dI.mb,dI.mr,dI.maa,dI.mbb,dI.mrr,dI.mab,dI.mar,dI.mbr,
                   dV.sa, dV.sb, dV.sr,
                   dT.s,
                   dE.sa,dE.sb,dE.sr,dE.saa,dE.sbb,dE.srr,dE.sab,dE.sar,dE.sbr,
                   dI.sa,dI.sb,dI.sr,dI.saa,dI.sbb,dI.srr,dI.sab,dI.sar,dI.sbr
                   ))
        })
    }
    out <- as.data.frame(ode(state,t,wv.BTV.coinfection.reassort.odes,parameters,
                             withlike=withlike))
    return(out)
}

rho.1 <- function(omega) {
    return(omega / 2 + (1 - omega) / 2 ^ 10)
}
rho.r <- function(omega) {
    return((1 - omega) * (1 - 1 / 2 ^ 9))
}

## As before, except allowing differing production rates between the viruses
wv.BTV.coinfection.reassort.variable <- function(t, state, parameters, withlike=0) {
    wv.BTV.coinfection.reassort.odes <- function(t, state, parameters, withlike=0)
    {
        with(as.list(c(state, parameters)),{
            ## Bloodmeal virus
            dV.ba = -beta.b * V.ba * (T.m + E.ma + E.mb + E.mr) - c.b * V.ba
            dV.bb = -beta.b * V.bb * (T.m + E.ma + E.mb + E.mr) - c.b * V.bb
            dV.br = -beta.b * V.br * (T.m + E.ma + E.mb + E.mr) - c.b * V.br
            
            ## Midgut virus
            dV.ma = (-beta.m * V.ma * (T.m + E.ma + E.mb + E.mr) + d.s * p.sa * (I.sa + I.saa) +
                     (1 - d.m) * p.ma * (I.ma + I.maa) - c.m * V.ma +
                     rho.1(withlike) * (d.s * p.sa * (I.sar + I.sab) +
                                        (1 - d.m) * p.ma * (I.mar + I.mab)))
            dV.mb = (-beta.m * V.mb * (T.m + E.ma + E.mb + E.mr) + d.s * p.sb* (I.sb + I.sbb) +
                     (1 - d.m) * p.mb * (I.mb + I.mbb) - c.m * V.mb +
                     rho.1(withlike) * (d.s * p.sb * (I.sbr + I.sab) +
                                        (1 - d.m) * p.mb * (I.mbr + I.mab)))
            dV.mr = (-beta.m * V.mr * (T.m + E.ma + E.mb + E.mr) + d.s * p.sr * (I.sr + I.srr) + 
                     (1 - d.m) * p.mr * (I.mr + I.mrr) - c.m * V.mr +
                     rho.1(withlike) * (d.s * p.sr * (I.sar + I.sbr) +
                                        (1 - d.m) * p.mr * (I.mar + I.mbr)) +
                     rho.r(withlike) * (d.s * p.sr * (I.sar + I.sbr + I.sab) +
                                        (1 - d.m) * p.mr * (I.mar + I.mbr + I.mab)))
            
            ## Midgut cells
            dT.m = -(beta.b * (V.ba + V.bb + V.br) +
                     beta.m * (V.ma + V.mb + V.mr)) * T.m + mu.m * (T0.m - T.m)
            
            dE.ma = (beta.b * V.ba + beta.m * V.ma) * T.m - (beta.b * V.bb + beta.m * V.mb) * E.ma -
                (beta.b * V.br + beta.m * V.mr) * E.ma - (beta.b * V.ba + beta.m * V.ma) * E.ma -
                E.ma / epsilon.m - mu.m * E.ma
            dE.mb = (beta.b * V.bb + beta.m * V.mb) * T.m - (beta.b * V.ba + beta.m * V.ma) * E.mb -
                (beta.b * V.br + beta.m * V.mr) * E.mb - (beta.b * V.bb + beta.m * V.mb) * E.mb -
                E.mb / epsilon.m - mu.m * E.mb
            dE.mr = (beta.b * V.br + beta.m * V.mr) * T.m - (beta.b * V.ba + beta.m * V.ma) * E.mr -
                (beta.b * V.bb + beta.m * V.mb) * E.mr - (beta.b * V.br + beta.m * V.mr) * E.mr -
                E.mr / epsilon.m - mu.m * E.mr
            dE.maa = (beta.b * V.ba + beta.m * V.ma) * E.ma - E.maa / epsilon.m - mu.m * E.maa
            dE.mbb = (beta.b * V.bb + beta.m * V.mb) * E.mb - E.mbb / epsilon.m - mu.m * E.mbb
            dE.mrr = (beta.b * V.br + beta.m * V.mr) * E.mr - E.mrr / epsilon.m - mu.m * E.mrr
            dE.mab = (beta.b * V.ba + beta.m * V.ma) * E.mb +
                (beta.b * V.bb + beta.m * V.mb) * E.ma - E.mab / epsilon.m - mu.m * E.mab
            dE.mar = (beta.b * V.ba + beta.m * V.ma) * E.mr +
                (beta.b * V.br + beta.m * V.mr) * E.ma - E.mar / epsilon.m - mu.m * E.mar
            dE.mbr = (beta.b * V.bb + beta.m * V.mb) * E.mr +
                (beta.b * V.br + beta.m * V.mr) * E.mb - E.mbr / epsilon.m - mu.m * E.mbr

            dI.ma = E.ma / epsilon.m - mui.m * I.ma
            dI.mb = E.mb / epsilon.m - mui.m * I.mb
            dI.mr = E.mr / epsilon.m - mui.m * I.mr
            dI.maa = E.maa / epsilon.m - mui.m * I.maa
            dI.mbb = E.mbb / epsilon.m - mui.m * I.mbb
            dI.mrr = E.mrr / epsilon.m - mui.m * I.mrr
            dI.mab = E.mab / epsilon.m - mui.m * I.mab
            dI.mar = E.mar / epsilon.m - mui.m * I.mar
            dI.mbr = E.mbr / epsilon.m - mui.m * I.mbr
            
            ## Secondary tissue virus
            dV.sa = (-beta.s * V.sa * (T.s + E.sa + E.sb + E.sr) + d.m * p.ma * (I.ma + I.maa) +
                (1 - d.s) * p.sa * (I.sa + I.saa) - c.s * V.sa +
                rho.1(withlike) * (d.m * p.ma * (I.mar + I.mab) + (1 - d.s) * p.sa * (I.sar + I.sab)))
            dV.sb = (-beta.s * V.sb * (T.s + E.sa + E.sb + E.sr) + d.m * p.mb * (I.mb + I.mbb) +
                (1 - d.s) * p.sb * (I.sb + I.sbb) - c.s * V.sb +
                rho.1(withlike) * (d.m * p.mb * (I.mbr + I.mab) + (1 - d.s) * p.sb * (I.sbr + I.sab)))
            dV.sr = (-beta.s * V.sr * (T.s + E.sa + E.sb + E.sr) + d.m * p.mr * (I.mr + I.mrr) +
                (1 - d.s) * p.sr * (I.sr + I.srr) - c.s * V.sr +
                rho.1(withlike) * (d.m * p.mr * (I.mar + I.mbr) + (1 - d.s) * p.sr * (I.sar + I.sbr)) +
                rho.r(withlike) * (d.m * p.mr * (I.mar + I.mbr + I.mab) +
                                   (1 - d.s) * p.sr * (I.sar + I.sbr + I.sab)))
            
            ## Secondary tissue cells
            dT.s = -beta.s * (V.sa + V.sb + V.sr) * T.s + mu.s * (T0.s - T.s)

            dE.sa = beta.s * V.sa * T.s - beta.s * V.sb * E.sa - beta.s * V.sr * E.sa -
                beta.s * V.sa * E.sa - E.sa / epsilon.s - mu.s * E.sa
            dE.sb = beta.s * V.sb * T.s - beta.s * V.sa * E.sb - beta.s * V.sr * E.sb -
                beta.s * V.sb * E.sb - E.sb / epsilon.s - mu.s * E.sb
            dE.sr = beta.s * V.sr * T.s - beta.s * V.sa * E.sr - beta.s * V.sb * E.sr -
                beta.s * V.sr * E.sr - E.sr / epsilon.s - mu.s * E.sr
            dE.saa = beta.s * V.sa * E.sa - E.saa / epsilon.s - mu.s * E.saa
            dE.sbb = beta.s * V.sb * E.sb - E.sbb / epsilon.s - mu.s * E.sbb
            dE.srr = beta.s * V.sr * E.sr - E.srr / epsilon.s - mu.s * E.srr
            dE.sab = beta.s * V.sa * E.sb + beta.s * V.sb * E.sa - E.sab / epsilon.s - mu.s * E.sab
            dE.sar = beta.s * V.sa * E.sr + beta.s * V.sr * E.sa - E.sar / epsilon.s - mu.s * E.sar
            dE.sbr = beta.s * V.sb * E.sr + beta.s * V.sr * E.sb - E.sbr / epsilon.s - mu.s * E.sbr
            
            dI.sa = E.sa / epsilon.s - mui.s * I.sa
            dI.sb = E.sb / epsilon.s - mui.s * I.sb
            dI.sr = E.sr / epsilon.s - mui.s * I.sr
            dI.saa = E.saa / epsilon.s - mui.s * I.saa
            dI.sbb = E.sbb / epsilon.s - mui.s * I.sbb
            dI.srr = E.srr / epsilon.s - mui.s * I.srr
            dI.sab = E.sab / epsilon.s - mui.s * I.sab
            dI.sar = E.sar / epsilon.s - mui.s * I.sar
            dI.sbr = E.sbr / epsilon.s - mui.s * I.sbr
            
            list(c(dV.ba, dV.bb, dV.br,
                   dV.ma, dV.mb, dV.mr,
                   dT.m,
                   dE.ma,dE.mb,dE.mr,dE.maa,dE.mbb,dE.mrr,dE.mab,dE.mar,dE.mbr,
                   dI.ma,dI.mb,dI.mr,dI.maa,dI.mbb,dI.mrr,dI.mab,dI.mar,dI.mbr,
                   dV.sa, dV.sb, dV.sr,
                   dT.s,
                   dE.sa,dE.sb,dE.sr,dE.saa,dE.sbb,dE.srr,dE.sab,dE.sar,dE.sbr,
                   dI.sa,dI.sb,dI.sr,dI.saa,dI.sbb,dI.srr,dI.sab,dI.sar,dI.sbr
                   ))
        })
    }
    out <- as.data.frame(ode(state,t,wv.BTV.coinfection.reassort.odes,parameters,
                             withlike=withlike))
    return(out)
}

## Reduced forms of the coinfection system, mainly for checking
## To obtain a model without reassortment, just set withlike = 1 in full model
wv.BTV.coinfection.simple <- function(t, state, parameters) {
    wv.BTV.coinfection.simple.odes <- function(t, state, parameters)
    {
        with(as.list(c(state, parameters)),{
            ## Bloodmeal virus
            dV.ba = -beta.b * V.ba * (T.m + E.ma + E.mb) - c.b * V.ba
            dV.bb = -beta.b * V.bb * (T.m + E.ma + E.mb) - c.b * V.bb
            
            ## Produced virus
            dV.ma = (-beta.m * V.ma * (T.m + E.mb + E.ma) +
                     p.m * (I.ma + I.maa) - c.m * V.ma +
                     rho.1(1) * p.m * I.mab)
            dV.mb = (-beta.m * V.mb * (T.m + E.ma + E.mb) +
                     p.m * (I.mb + I.mbb) - c.m * V.mb +
                     rho.1(1) * p.m * I.mab)
            
            ## Cells
            dT.m = -(beta.b * (V.ba + V.bb) +
                     beta.m * (V.ma + V.mb)) * T.m + mu.m * (T0.m - T.m)
            
            dE.ma = (beta.b * V.ba + beta.m * V.ma) * T.m - (beta.b * V.bb + beta.m * V.mb) * E.ma -
                (beta.b * V.ba + beta.m * V.ma) * E.ma - E.ma / epsilon.m - mu.m * E.ma
            dE.mb = (beta.b * V.bb + beta.m * V.mb) * T.m - (beta.b * V.ba + beta.m * V.ma) * E.mb -
                (beta.b * V.bb + beta.m * V.mb) * E.mb - E.mb / epsilon.m - mu.m * E.mb
            dE.maa = (beta.b * V.ba + beta.m * V.ma) * E.ma - E.maa / epsilon.m - mu.m * E.maa
            dE.mbb = (beta.b * V.bb + beta.m * V.mb) * E.mb - E.mbb / epsilon.m - mu.m * E.mbb
            dE.mab = (beta.b * V.ba + beta.m * V.ma) * E.mb +
                (beta.b * V.bb + beta.m * V.mb) * E.ma - E.mab / epsilon.m - mu.m * E.mab

            dI.ma = E.ma / epsilon.m - mui.m * I.ma
            dI.mb = E.mb / epsilon.m - mui.m * I.mb
            dI.maa = E.maa / epsilon.m - mui.m * I.maa
            dI.mbb = E.mbb / epsilon.m - mui.m * I.mbb
            dI.mab = E.mab / epsilon.m - mui.m * I.mab
                        
            list(c(dV.ba, dV.bb,
                   dV.ma, dV.mb,
                   dT.m,
                   dE.ma,dE.mb,dE.maa,dE.mbb,dE.mab,
                   dI.ma,dI.mb,dI.maa,dI.mbb,dI.mab))
        })
    }
    out <- as.data.frame(ode(state,t,wv.BTV.coinfection.simple.odes,parameters))
    return(out)
}

wv.BTV.coinfection.simple.nonneutral <- function(t, state, parameters) {
    wv.BTV.coinfection.simple.odes <- function(t, state, parameters)
    {
        with(as.list(c(state, parameters)),{
            ## Bloodmeal virus
            dV.ba = -beta.b * V.ba * (T.m + E.mb) - c.b * V.ba
            dV.bb = -beta.b * V.bb * (T.m + E.ma) - c.b * V.bb
            
            ## Produced virus
            dV.ma = (-beta.m * V.ma * (T.m + E.mb) +
                     p.m * I.ma - c.m * V.ma +
                     rho.1(1) * p.m * I.mab)
            dV.mb = (-beta.m * V.mb * (T.m + E.ma) +
                     p.m * I.mb - c.m * V.mb +
                     rho.1(1) * p.m * I.mab)
            
            ## Cells
            dT.m = -(beta.b * (V.ba + V.bb) +
                     beta.m * (V.ma + V.mb)) * T.m + mu.m * (T0.m - T.m)
            
            dE.ma = (beta.b * V.ba + beta.m * V.ma) * T.m - (beta.b * V.bb + beta.m * V.mb) * E.ma -
                E.ma / epsilon.m - mu.m * E.ma
            dE.mb = (beta.b * V.bb + beta.m * V.mb) * T.m - (beta.b * V.ba + beta.m * V.ma) * E.mb -
                E.mb / epsilon.m - mu.m * E.mb
            dE.mab = (beta.b * V.ba + beta.m * V.ma) * E.mb +
                (beta.b * V.bb + beta.m * V.mb) * E.ma - E.mab / epsilon.m - mu.m * E.mab

            dI.ma = E.ma / epsilon.m - mui.m * I.ma
            dI.mb = E.mb / epsilon.m - mui.m * I.mb
            dI.mab = E.mab / epsilon.m - mui.m * I.mab
                        
            list(c(dV.ba, dV.bb,
                   dV.ma, dV.mb,
                   dT.m,
                   dE.ma,dE.mb,dE.mab,
                   dI.ma,dI.mb,dI.mab))

        })
    }
    out <- as.data.frame(ode(state,t,wv.BTV.coinfection.simple.odes,parameters))
    return(out)
}

## The calculation of the output in the coinfection model with reassortment
calc.R.full <- function(state.coinf,parms,fitting.parms) {
    ## Calculate each spacing between infections separately.
    ## When spacing != 0, decrease the p.db by a fixed amount
    ## Do this on a logistic scale so that it doesn't go above 1
    ## But ensure that the probability increase only takes place subsequent to second bloodmeal
    withlike <- fitting.parms["withlike"]
    ##parms["epsilon.m"] <- fitting.parms["epsilon.m"]
    ##parms["epsilon.s"] <- fitting.parms["epsilon.s"]
    output.list <- list()
    for (second.intro in unique(combined.runtimes$gap)){
        runtime <- combined.runtimes[gap==second.intro,runtime]
        if (second.intro > 0) {
            p.db.new <- plogis(qlogis(p.db) + fitting.parms[["p.db.inc"]])
            times.a <- seq(0,second.intro,1)
            times.b <- seq(second.intro,runtime,1)
            ## 1. no infection 1 - p.mib
            state.1a <- state.coinf
            state.1a["T.m"] <- 0
            state.1a["T.s"] <- 0
            state.1a["V.ba"] <- floor(initial.titre*log(2) + 0.5)
            parms.1 <- parms
            parms.1["T0.m"] <- 0
            parms.1["T0.s"] <- 0
            dat.1a = wv.BTV.coinfection.reassort(times.a,state.1a,parms.1,withlike)
            state.1b <- pmax(as.numeric(dat.1a[nrow(dat.1a),-1]),0)
            names(state.1b) <- names(state.1a)
            state.1b["V.bb"] <- floor(initial.titre*log(2) + 0.5)
            dat.1b = wv.BTV.coinfection.reassort(times.b,state.1b,parms.1,withlike)
            dat.1 <- rbind(dat.1a,dat.1b[-1,])
            V.tot.1 <- rowSums(dat.1[,grepl("V.",names(dat.1))])
            R.tot.1 <- rowSums(dat.1[,c("V.br","V.mr","V.sr")])

            ## 2. constrained to midgut p.mib(1 - p.db.new)
            state.2a <- state.1a
            state.2a["T.m"] <- parms[["T0.m"]]
            parms.2 <- parms
            parms.2["T0.s"] <- 0
            dat.2a <- wv.BTV.coinfection.reassort(times.a,state.2a,parms.2,withlike)
            state.2b <- pmax(as.numeric(dat.2a[nrow(dat.2a),-1]),0)
            names(state.2b) <- names(state.2a)
            state.2b[["V.bb"]] <- floor(initial.titre*log(2) + 0.5)
            dat.2b <- wv.BTV.coinfection.reassort(times.b,state.2b,parms.2,withlike)
            dat.2 <- rbind(dat.2a,dat.2b[-1,])
            V.tot.2 <- rowSums(dat.2[,grepl("V.",names(dat.2))])
            R.tot.2 <- rowSums(dat.2[,c("V.br","V.mr","V.sr")])

            ## 3. disseminated infection following second bloodmeal p.mib*(p.db.new - p.db) 
            state.3a <- state.2a
            parms.3a <- parms
            parms.3a["T0.s"] <- 0
            dat.3a <- wv.BTV.coinfection.reassort(times.a,state.3a,parms.3a,withlike)
            state.3b <- pmax(as.numeric(dat.3a[nrow(dat.3a),-1]),0)
            names(state.3b) <- names(state.3a)
            state.3b[["V.bb"]] <- floor(initial.titre*log(2) + 0.5)
            state.3b["T.s"] <- parms[["T0.s"]]
            parms.3b <- parms
            dat.3b <- wv.BTV.coinfection.reassort(times.b,state.3b,parms.3b,withlike)
            dat.3 <- rbind(dat.3a,dat.3b[-1,])
            V.tot.3 <- rowSums(dat.3[,grepl("V.",names(dat.3))])
            R.tot.3 <- rowSums(dat.3[,c("V.br","V.mr","V.sr")])

            ## 4. disseminated infection right away p.mib*p.db 
            state.4a <- state.3a
            state.4a["T.s"] <- parms[["T0.s"]]
            parms.4 <- parms
            dat.4a <- wv.BTV.coinfection.reassort(times.a,state.4a,parms.4,withlike)
            state.4b <- pmax(as.numeric(dat.4a[nrow(dat.4a),-1]),0)
            names(state.4b) <- names(state.4a)
            state.4b[["V.bb"]] <- floor(initial.titre*log(2) + 0.5)
            dat.4b <- wv.BTV.coinfection.reassort(times.b,state.4b,parms.4,withlike)
            dat.4 <- rbind(dat.4a,dat.4b[-1,])
            V.tot.4 <- rowSums(dat.4[,grepl("V.",names(dat.4))])
            R.tot.4 <- rowSums(dat.4[,c("V.br","V.mr","V.sr")])

            ## Get weighted sum
            V.tot <- V.tot.1*(1-p.mib) + V.tot.2*p.mib*(1-p.db.new) +
                V.tot.3*p.mib*(p.db.new-p.db) + V.tot.4*p.mib*p.db
            V.tot.pos <- V.tot.2*(1-p.db.new) + V.tot.3*(p.db.new-p.db) + V.tot.4*p.db
            R.tot <- R.tot.1*(1-p.mib) + R.tot.2*p.mib*(1-p.db.new) +
                R.tot.3*p.mib*(p.db.new-p.db) + R.tot.4*p.mib*p.db
            R.tot.pos <- R.tot.2*(1-p.db.new) + R.tot.3*(p.db.new-p.db) + R.tot.4*p.db
            output.list[[as.character(second.intro)]] <- list(V.tot.1=V.tot.1,
                                                              V.tot.2=V.tot.2,V.tot.3=V.tot.3,
                                                              V.tot.4=V.tot.4,
                                                              V.tot=V.tot,V.tot.pos=V.tot.pos,
                                                              R.tot.1=R.tot.1,R.tot.2=R.tot.2,
                                                              R.tot.3=R.tot.3,
                                                              R.tot.4=R.tot.4,
                                                              R.tot=R.tot,R.tot.pos=R.tot.pos)
        } else {
            times <- seq(0,runtime,1)
            ## 1. no infection
            state.1 <- state.coinf
            state.1[["V.ba"]] <- 0.5 * floor(initial.titre*log(2) + 0.5)
            state.1[["V.bb"]] <- 0.5 * floor(initial.titre*log(2) + 0.5)
            state.1["T.m"] <- 0
            state.1["T.s"] <- 0
            parms.1 <- parms
            parms.1["T0.m"] <- 0
            parms.1["T0.s"] <- 0
            dat.1 <- wv.BTV.coinfection.reassort(times,state.1,parms.1,withlike)
            V.tot.1 <- rowSums(dat.1[,grepl("V.",names(dat.1))])
            R.tot.1 <- rowSums(dat.1[,c("V.br","V.mr","V.sr")])

            ## 2. constrained to midgut p.mib(1 - p.db)
            state.2 <- state.1
            state.2["T.m"] <- parms[["T0.m"]]
            parms.2 <- parms
            parms.2["T0.s"] <- 0
            dat.2 <- wv.BTV.coinfection.reassort(times,state.2,parms.2,withlike)
            V.tot.2 <- rowSums(dat.2[,grepl("V.",names(dat.2))])
            R.tot.2 <- rowSums(dat.2[,c("V.br","V.mr","V.sr")])

            ## 3. disseminated infection p.mib*p.db
            state.3 <- state.2
            state.3["T.s"] <- parms[["T0.s"]]
            parms.3 <- parms
            dat.3 <- wv.BTV.coinfection.reassort(times,state.3,parms.3,withlike)
            V.tot.3 <- rowSums(dat.3[,grepl("V.",names(dat.3))])
            R.tot.3 <- rowSums(dat.3[,c("V.br","V.mr","V.sr")])

            ## Get weighted sum
            V.tot <- V.tot.1*(1-p.mib) + V.tot.2*p.mib*(1-p.db) +
                V.tot.3*p.mib*p.db
            V.tot.pos <- V.tot.2*(1-p.db) + V.tot.3*p.db
            R.tot <- R.tot.1*(1-p.mib) + R.tot.2*p.mib*(1-p.db) +
                R.tot.3*p.mib*p.db
            R.tot.pos <- R.tot.2*(1-p.db) + R.tot.3*p.db
            output.list[[as.character(second.intro)]] <- list(V.tot.1=V.tot.1,V.tot.2=V.tot.2,
                                                              V.tot.3=V.tot.3,
                                                              V.tot=V.tot,V.tot.pos=V.tot.pos,
                                                              R.tot.1=R.tot.1,R.tot.2=R.tot.2,
                                                              R.tot.3=R.tot.3,
                                                              R.tot=R.tot,R.tot.pos=R.tot.pos)
        }
    }
    return(output.list)
}

calc.R.variable <- function(state.coinf,parms,fitting.parms,runtime,second.intro) {
    ## Calculate each spacing between infections separately.
    ## When spacing != 0, decrease the p.db by a fixed amount
    ## Do this on a logistic scale so that it doesn't go above 1
    ## But ensure that the probability increase only takes place subsequent to second bloodmeal
    withlike <- fitting.parms["withlike"]
    ##parms["epsilon.m"] <- fitting.parms["epsilon.m"]
    ##parms["epsilon.s"] <- fitting.parms["epsilon.s"]
    if (second.intro > 0) {
        p.db.new <- plogis(qlogis(p.db) + fitting.parms[["p.db.inc"]])
        times.a <- seq(0,second.intro,1)
        times.b <- seq(second.intro,runtime,1)
        ## 1. no infection 1 - p.mib
        state.1a <- state.coinf
        state.1a["T.m"] <- 0
        state.1a["T.s"] <- 0
        parms.1 <- parms
        parms.1["T0.m"] <- 0
        parms.1["T0.s"] <- 0
        dat.1a = wv.BTV.coinfection.reassort.variable(times.a,state.1a,parms.1,withlike)
        state.1b <- pmax(as.numeric(dat.1a[nrow(dat.1a),-1]),0)
        names(state.1b) <- names(state.1a)
        state.1b["V.bb"] <- floor(initial.titre*log(2) + 0.5)
        dat.1b = wv.BTV.coinfection.reassort.variable(times.b,state.1b,parms.1,withlike)
        dat.1 <- rbind(dat.1a,dat.1b[-1,])
        V.tot.1 <- rowSums(dat.1[,grepl("V.",names(dat.1))])
        R.tot.1 <- rowSums(dat.1[,c("V.br","V.mr","V.sr")])

        ## 2. constrained to midgut p.mib(1 - p.db.new)
        state.2a <- state.1a
        state.2a["T.m"] <- parms[["T0.m"]]
        parms.2 <- parms
        parms.2["T0.s"] <- 0
        dat.2a <- wv.BTV.coinfection.reassort.variable(times.a,state.2a,parms.2,withlike)
        state.2b <- pmax(as.numeric(dat.2a[nrow(dat.2a),-1]),0)
        names(state.2b) <- names(state.2a)
        state.2b[["V.bb"]] <- floor(initial.titre*log(2) + 0.5)
        dat.2b <- wv.BTV.coinfection.reassort.variable(times.b,state.2b,parms.2,withlike)
        dat.2 <- rbind(dat.2a,dat.2b[-1,])
        V.tot.2 <- rowSums(dat.2[,grepl("V.",names(dat.2))])
        R.tot.2 <- rowSums(dat.2[,c("V.br","V.mr","V.sr")])

        ## 3. disseminated infection following second bloodmeal p.mib*(p.db.new - p.db) 
        state.3a <- state.2a
        parms.3a <- parms
        parms.3a["T0.s"] <- 0
        dat.3a <- wv.BTV.coinfection.reassort.variable(times.a,state.3a,parms.3a,withlike)
        state.3b <- pmax(as.numeric(dat.3a[nrow(dat.3a),-1]),0)
        names(state.3b) <- names(state.3a)
        state.3b[["V.bb"]] <- floor(initial.titre*log(2) + 0.5)
        state.3b["T.s"] <- parms[["T0.s"]]
        parms.3b <- parms
        dat.3b <- wv.BTV.coinfection.reassort.variable(times.b,state.3b,parms.3b,withlike)
        dat.3 <- rbind(dat.3a,dat.3b[-1,])
        V.tot.3 <- rowSums(dat.3[,grepl("V.",names(dat.3))])
        R.tot.3 <- rowSums(dat.3[,c("V.br","V.mr","V.sr")])

        ## 4. disseminated infection right away p.mib*p.db 
        state.4a <- state.3a
        state.4a["T.s"] <- parms[["T0.s"]]
        parms.4 <- parms
        dat.4a <- wv.BTV.coinfection.reassort.variable(times.a,state.4a,parms.4,withlike)
        state.4b <- pmax(as.numeric(dat.4a[nrow(dat.4a),-1]),0)
        names(state.4b) <- names(state.4a)
        state.4b[["V.bb"]] <- floor(initial.titre*log(2) + 0.5)
        dat.4b <- wv.BTV.coinfection.reassort.variable(times.b,state.4b,parms.4,withlike)
        dat.4 <- rbind(dat.4a,dat.4b[-1,])
        V.tot.4 <- rowSums(dat.4[,grepl("V.",names(dat.4))])
        R.tot.4 <- rowSums(dat.4[,c("V.br","V.mr","V.sr")])

        ## Get weighted sum
        V.tot <- V.tot.1*(1-p.mib) + V.tot.2*p.mib*(1-p.db.new) +
            V.tot.3*p.mib*(p.db.new-p.db) + V.tot.4*p.mib*p.db
        V.tot.pos <- V.tot.2*(1-p.db.new) + V.tot.3*(p.db.new-p.db) + V.tot.4*p.db
        R.tot <- R.tot.1*(1-p.mib) + R.tot.2*p.mib*(1-p.db.new) +
            R.tot.3*p.mib*(p.db.new-p.db) + R.tot.4*p.mib*p.db
        R.tot.pos <- R.tot.2*(1-p.db.new) + R.tot.3*(p.db.new-p.db) + R.tot.4*p.db
        output.list <- list(V.tot=V.tot,V.tot.pos=V.tot.pos,
                            R.tot=R.tot,R.tot.pos=R.tot.pos)
    } else {
        times <- seq(0,runtime,1)
        ## 1. no infection
        state.1 <- state.coinf
        state.1[["V.ba"]] <- 0.5 * floor(initial.titre*log(2) + 0.5)
        state.1[["V.bb"]] <- 0.5 * floor(initial.titre*log(2) + 0.5)
        state.1["T.m"] <- 0
        state.1["T.s"] <- 0
        parms.1 <- parms
        parms.1["T0.m"] <- 0
        parms.1["T0.s"] <- 0
        dat.1 <- wv.BTV.coinfection.reassort.variable(times,state.1,parms.1,withlike)
        V.tot.1 <- rowSums(dat.1[,grepl("V.",names(dat.1))])
        R.tot.1 <- rowSums(dat.1[,c("V.br","V.mr","V.sr")])

        ## 2. constrained to midgut p.mib(1 - p.db)
        state.2 <- state.1
        state.2["T.m"] <- parms[["T0.m"]]
        parms.2 <- parms
        parms.2["T0.s"] <- 0
        dat.2 <- wv.BTV.coinfection.reassort.variable(times,state.2,parms.2,withlike)
        V.tot.2 <- rowSums(dat.2[,grepl("V.",names(dat.2))])
        R.tot.2 <- rowSums(dat.2[,c("V.br","V.mr","V.sr")])

        ## 3. disseminated infection p.mib*p.db
        state.3 <- state.2
        state.3["T.s"] <- parms[["T0.s"]]
        parms.3 <- parms
        dat.3 <- wv.BTV.coinfection.reassort.variable(times,state.3,parms.3,withlike)
        V.tot.3 <- rowSums(dat.3[,grepl("V.",names(dat.3))])
        R.tot.3 <- rowSums(dat.3[,c("V.br","V.mr","V.sr")])

        ## Get weighted sum
        V.tot <- V.tot.1*(1-p.mib) + V.tot.2*p.mib*(1-p.db) +
            V.tot.3*p.mib*p.db
        V.tot.pos <- V.tot.2*(1-p.db) + V.tot.3*p.db
        R.tot <- R.tot.1*(1-p.mib) + R.tot.2*p.mib*(1-p.db) +
            R.tot.3*p.mib*p.db
        R.tot.pos <- R.tot.2*(1-p.db) + R.tot.3*p.db
        output.list <- list(V.tot=V.tot,V.tot.pos=V.tot.pos,
                            R.tot=R.tot,R.tot.pos=R.tot.pos)
    }    
    return(output.list)
}

## The within-vector model in Tuncer et al.
## https://www.tandfonline.com/doi/full/10.1080/17513758.2021.1970261
tuncer.wv = function(t, state, parameters) 
{
    wv.odes <- function(t,state,parameters) {
        with(as.list(c(state, parameters)),{
            dPm = rm * Pm * (1 - Pm / Km) - mt * Pm
            dPs = mt * Pm^2 / (B + Pm^2) + rs * Ps * (1 - Ps /Ks)
            list(c(dPm, dPs))
        })
    }
    out <- as.data.frame(ode(state,t,wv.odes,parameters))
    return(out)
}

## parms <- c(rm=1.26,Km=4.33e5,mt=0.11,B=1.95e4,rs=1.77,Ks=6.84e6)
## state <- c(Pm=1e3,Ps=0)
