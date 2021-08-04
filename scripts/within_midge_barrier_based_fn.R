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
wv.BTV.barrier.det = function(t, state, parameters) {
    wv.BTV.odes <- function(t, state, parameters)
        with(as.list(c(state, parameters)),{
            dV.m = -c.m * V.m - beta.m * V.m * T.m
            dT.m = -beta.m * V.m * T.m 
            dE.m = beta.m * V.m * T.m - E.m / epsilon.m
            dI.m = E.m / epsilon.m
            dV.s = p.m * I.m + p.s * I.s - c.s * V.s - beta.s * T.s * V.s
            dT.s = -beta.s * T.s * V.s 
            dE.s = beta.s * T.s * V.s - E.s / epsilon.s
            dI.s = E.s / epsilon.s
            list(c(dV.m, dT.m, dE.m, dI.m, 
                   dV.s, dT.s, dE.s, dI.s))
        })

    out <- as.data.frame(ode(state,t,wv.BTV.odes,parameters))
    return(out)
}

## stochastic model
wv.BTV.barrier.stoch <- function(t, state, parameters){

    ## Transitions matrix
    wv.BTV.transitions <- list(
        c(V.m=-1),
        c(T.m=-1, V.m=-1, E.m=1),
        c(E.m=-1,I.m=1),
        c(V.s=1),
        c(V.s=-1),
        c(T.s=-1, V.s=-1, E.s=1),
        c(E.s=-1,I.s=1)
    )

    wv.BTV.rates <- function(state, parameters, t)
    {
        with(as.list(c(state, parameters)),{
            return(c(c.m * V.m,
                     beta.m * V.m * T.m,
                     E.m / epsilon.m,
                     p.m * I.m + p.s * I.s,
                     c.s * V.s,
                     beta.s * T.s * V.s,
                     E.s/epsilon.s
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

wv.BTV.coinfection.reassort = function(t, state, parameters, prod.fn.a=production.fn.a,
                                       prod.fn.b=production.fn.b, prod.fn.r=production.fn.r)
{
    with(as.list(c(state, parameters)),{
        ## Stuff to adapt:
        
        dT.m = -(beta.b * V.b + beta.m * V.m) * T.m 
        dE.m1 = (beta.b * V.b + beta.m * V.m) * T.m - 2 * E.m1 / epsilon.m
        dE.m2 = 2 * (E.m1 - E.m2) / epsilon.m 
        dI.m = 2 * E.m2 / epsilon.m
        dV.l = d.m * p.m * I.m + (1 - d.l) * p.l * I.l - c.l * V.l
        dT.l = -beta.l * T.l * V.l 
        dE.l1 = beta.l * T.l * V.l  - 2 * E.l1 / epsilon.l
        dE.l2 = 2 * (E.l1 - E.l2) / epsilon.l 
        dI.l = 2 * E.l2 / epsilon.l
        dV.d = d.l * p.l * I.l + (1 - d.d) * p.d * I.d - c.d * V.d
        dT.d = -beta.d * T.d * V.d 
        dE.d1 = beta.d * T.d * V.d - 2 * E.d1 / epsilon.d
        dE.d2 = 2 * (E.d1 - E.d2) / epsilon.d
        dI.d = 2 * E.d2 / epsilon.d
        dV.s = d.d * p.d * I.d + p.s * I.s - c.s * V.s
        dT.s = -beta.s * T.s * V.s 
        dE.s1 = beta.s * T.s * V.s - 2 * E.s1 / epsilon.s
        dE.s2 = 2 * (E.s1 - E.s2) / epsilon.s
        dI.s = 2 * E.s2 / epsilon.s

        ## Bloodmeal virus
        dV.ba = -c.m * V.ba
        dV.bb = -c.m * V.bb
        dV.br = -c.m * V.br
        ## Midgut virus START HERE
        dV.ma = (1 - d.m) * p.m * I.m - c.m * V.m
        dV.mb = 
        dV.mr =
        ## Midgut cells
        dT.m = -beta.m * T.m * (V.ma + V.mb)
        dI.ma0 = beta.m * T.m * V.ma - I.ma0 * (n / epsilon + beta.m * V.mb)
        dI.mb0 = beta.m * T.m * V.mb - I.mb0 * (n / epsilon + beta.m * V.ma)
        dI.man = n * get(paste0("I.ma",n-1)) / epsilon
        dI.mbn = n * get(paste0("I.mb",n-1)) / epsilon
        dI.ma <- vector(mode="numeric",length=n-1)
        dI.mb <- vector(mode="numeric",length=n-1)
        dI.mabj0 = vector(mode="numeric",length=n-1)
        dI.mab0k = vector(mode="numeric",length=n-1)
        dI.mabjn = vector(mode="numeric",length=n-1)
        dI.mabnk = vector(mode="numeric",length=n-1)
        ## Secondary tissue virus
        sumV.ma.vec <- vector(mode="numeric",length=n-1)
        sumV.mb.vec <- vector(mode="numeric",length=n-1)
        sumV.mr.vec <- vector(mode="numeric",length=n-1)
        ## Secondary tissue cells
        dT.s = -beta.s * T.s * (V.sa + V.sb + V.sr)
        dI.sa0 = beta.s * T.s * V.sa - I.sa0 * (1 / epsilon + beta.s * (V.sb + V.sr))
        dI.sb0 = beta.s * T.s * V.sb - I.sb0 * (1 / epsilon + beta.s * (V.sa + V.sr))
        dI.sr0 = beta.s * T.s * V.sr - I.sr0 * (1 / epsilon + beta.s * (V.sa + V.sb))
        dI.san = I.sa0 / epsilon
        dI.sbn = I.sb0 / epsilon
        dI.srn = I.sr0 / epsilon
        dI.sab00 <- beta.s * (V.sa * I.sb0 + V.sb * I.sa0) - I.sab00/epsilon
        dI.sar00 <- beta.s * (V.sa * I.sr0 + V.sr * I.sa0) - I.sar00/epsilon
        dI.sbr00 <- beta.s * (V.sr * I.sb0 + V.sb * I.sr0) - I.sbr00/epsilon
        dI.sabnn <- I.sab00/epsilon
        dI.sarnn <- I.sar00/epsilon
        dI.sbrnn <- I.sbr00/epsilon
        if (n>1) {
            for (ii in 1:(n-1)) {
                ## Midgut cells
                diffval.a <- get(paste0("I.ma",ii-1)) - get(paste0("I.ma",ii))
                diffval.b <- get(paste0("I.mb",ii-1)) - get(paste0("I.mb",ii))
                dI.ma[ii] <- n * diffval.a / epsilon - beta.m * V.mb * get(paste0("I.ma",ii))
                dI.mb[ii] <- n * diffval.b / epsilon - beta.m * V.ma * get(paste0("I.mb",ii))
                dI.mabj0[ii] <- beta.m * V.mb * get(paste0("I.ma",ii)) - n * get(paste0("I.mab",ii,"0")) / epsilon
                dI.mab0k[ii] <- beta.m * V.ma * get(paste0("I.mb",ii)) - n * get(paste0("I.mab0",ii)) / epsilon
                dI.mabjn[ii] <- n * get(paste0("I.mab",ii-1,n-1)) / epsilon
                dI.mabnk[ii] <- n * get(paste0("I.mab",n-1,ii-1)) / epsilon
                ## Secondary tissue virus
                sumV.ma.vec[ii] <- get(paste0("I.mab",ii,"n")) * prod.fn.a(ii,n,n) + get(paste0("I.mabn",ii)) * prod.fn.a(n,ii,n)
                sumV.mb.vec[ii] <- get(paste0("I.mab",ii,"n")) * prod.fn.b(ii,n,n) + get(paste0("I.mabn",ii)) * prod.fn.b(n,ii,n)
                sumV.mr.vec[ii] <- (get(paste0("I.mab",ii,"n")) + get(paste0("I.mabn",ii))) * prod.fn.r(ii,n,n)
            }
        }
        ## Midgut cells
        dI.mabjk = array(dim = c(n-1,n-1))
        if (n>1){
            for (jj in 1:(n-1)) {
                for (kk in 1:(n-1)) {
                    ## Midgut cells
                    diffval <- get(paste0("I.mab",jj-1,kk-1)) - get(paste0("I.mab",jj,kk))
                    dI.mabjk[jj,kk] <- n * diffval / epsilon
                }
            }
        }
        ## Midgut cells
        dI.mab00 <- beta.m * (V.ma * I.mb0 + V.mb * I.ma0) - n * I.mab00 / epsilon
        dI.mabnn <- n * get(paste0("I.mab",n-1,n-1)) / epsilon                
        ## Secondary tissue virus
        sumV.ma <- sum(sumV.ma.vec)
        sumV.mb <- sum(sumV.mb.vec)
        sumV.mr <- sum(sumV.mr.vec)
        ## Next three lines need editting to account for n=1
        if (n>1) {
            dV.sa = p.m * (I.man + sumV.ma + I.mabnn * prod.fn.a(n,n,n)) + k * p.s * (I.san + (I.sabnn + I.sarnn) * prod.fn.a(1,1,1)) - c.s * V.sa
            dV.sb = p.m * (I.mbn + sumV.mb + I.mabnn * prod.fn.b(n,n,n)) + k * p.s * (I.sbn + I.sabnn * prod.fn.b(1,1,1) + I.sbrnn * prod.fn.a(1,1,1)) - c.s * V.sb
            dV.sr = p.m * (sumV.mr + I.mabnn * prod.fn.r(n,n,n)) +
                k * p.s * ((I.sabnn + I.sarnn + I.sbrnn) * prod.fn.r(1,1,1) + I.srn +
                          (I.sarnn + I.sbrnn) * prod.fn.b(1,1,1))  - c.s * V.sr
        } else {
            dV.sa = p.m * (I.man + I.mabnn * prod.fn.a(n,n,n)) + k * p.s * (I.san + (I.sabnn + I.sarnn) * prod.fn.a(1,1,1)) - c.s * V.sa
            dV.sb = p.m * (I.mbn + I.mabnn * prod.fn.b(n,n,n)) + k * p.s * (I.sbn + I.sabnn * prod.fn.b(1,1,1) + I.sbrnn * prod.fn.a(1,1,1)) - c.s * V.sb
            dV.sr = p.m * I.mabnn * prod.fn.r(n,n,n) +
                k * p.s * ((I.sabnn + I.sarnn + I.sbrnn) * prod.fn.r(1,1,1) + I.srn +
                          (I.sarnn + I.sbrnn) * prod.fn.b(1,1,1))  - c.s * V.sr
        }
        list(c(dV.ma, dV.mb,
               dT.m,
               dI.ma0, dI.ma, dI.man,
               dI.mb0, dI.mb, dI.mbn,
               dI.mab00, dI.mabj0, dI.mab0k,
               as.numeric(dI.mabjk),
               dI.mabnk, dI.mabjn, dI.mabnn,
               dV.sa, dV.sb, dV.sr,
               dT.s,
               dI.sa0, dI.san,
               dI.sb0, dI.sbn,
               dI.sr0, dI.srn,
               dI.sab00,
               dI.sar00,
               dI.sbr00,
               dI.sabnn,
               dI.sarnn,
               dI.sbrnn))
    })
}
