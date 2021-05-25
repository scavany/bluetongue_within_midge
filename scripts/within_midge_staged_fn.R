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
    dI.s0 = beta.s * T.s * V.h - I.s0 / epsilon
    dI.s <- I.s0 / epsilon
    list(c(dV.m, dT.m, dI.m0, dI.m, dI.mn,
           dV.h, dT.s, dI.s0,dI.s))
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

## within midge model variable temperature NEEDS EDITING
wv.BTV.varT = function(t, state, parameters, incubation.period.fn)
{
    with(as.list(c(state, parameters)),{
        eip <- incubation.period.fn(t)
        tau <- eip - epsilon
        dV.m = -c.m * V.m
        dT.m = -beta.m * T.m * V.m
        dI.m0 = beta.m * T.m * V.m - n * I.m0 / epsilon
        dI.m <- vector(mode="numeric",length=n-1)
        for (i in 1:(n-1)) {
            diffval <- get(paste0("I.m",i-1)) - get(paste0("I.m",i))
            dI.m[i] <- n * diffval / epsilon
        }
        dI.mn = n * get(paste0("I.m",n-1)) / epsilon
        ## dV.h = p * I.mn - c.s * V.h # only midgut infected cells produce virus
        dV.h = p.m * I.mn + p.s * I.sn - c.s * V.h # all infected cells produce virus
        dT.s = -beta.s * T.s * V.h
        dI.s0 = beta.s * T.s * V.h - n * I.s0 / tau
        dI.s <- vector(mode="numeric",length=n-1)
        for (i in 1:(n-1)) {
            diffval <- get(paste0("I.s",i-1)) - get(paste0("I.s",i))
            dI.s[i] <- n * diffval / tau
        }
        dI.sn = n * get(paste0("I.s",n-1)) / tau
        list(c(dV.m, dT.m, dI.m0, dI.m, dI.mn,
               dV.h, dT.s, dI.s0, dI.s, dI.sn))
    })
}

## within midge reassortment
wv.BTV.coinfection.reassort = function(t, state, parameters, prod.fn.a=production.fn.a,
                                       prod.fn.b=production.fn.b, prod.fn.r=production.fn.r)
{
    with(as.list(c(state, parameters)),{
        ## Midgut virus
        dV.ma = -c.m * V.ma
        dV.mb = -c.m * V.mb
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
        if (n>1){
            dI.mabjk = array(dim = c(n-1,n-1))
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
               dI.sbr00))
    })
}

## NEEDS EDITING
wv.BTV.coinfection.reassort.varT = function(t, state, parameters, incubation.period.fn,
                                            prod.fn.a=production.fn.a,
                                            prod.fn.b=production.fn.b, prod.fn.r=production.fn.r)
{
    with(as.list(c(state, parameters)),{
        eip <- incubation.period.fn(t)
        tau <- eip - epsilon 
        ## Midgut virus
        dV.ma = -c.m * V.ma
        dV.mb = -c.m * V.mb
        ## Midgut cells
        dT.m = -beta.m * T.m * (V.ma + V.mb)
        dI.ma0 = beta.m * T.m * V.ma - I.ma0 * (n / epsilon + beta * V.mb)
        dI.mb0 = beta.m * T.m * V.mb - I.mb0 * (n / epsilon + beta * V.ma)
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
        sumV.sa.vec <- vector(mode="numeric",length=n-1)
        sumV.sb.vec <- vector(mode="numeric",length=n-1)
        sumV.sr.vec <- vector(mode="numeric",length=n-1)
        ## Secondary tissue cells
        dT.s = -beta.s * T.s * (V.sa + V.sb + V.sr)
        dI.sa0 = beta.s * T.s * V.sa - I.sa0 * (n / tau + beta * (V.sb + V.sr))
        dI.sb0 = beta.s * T.s * V.sb - I.sb0 * (n / tau + beta * (V.sa + V.sr))
        dI.sr0 = beta.s * T.s * V.sr - I.sr0 * (n / tau + beta * (V.sa + V.sb))
        dI.san = n * get(paste0("I.sa",n-1)) / tau
        dI.sbn = n * get(paste0("I.sb",n-1)) / tau
        dI.srn = n * get(paste0("I.sr",n-1)) / tau
        dI.sa <- vector(mode="numeric",length=n-1)
        dI.sb <- vector(mode="numeric",length=n-1)
        dI.sr <- vector(mode="numeric",length=n-1)
        dI.sabj0 = vector(mode="numeric",length=n-1)
        dI.sab0k = vector(mode="numeric",length=n-1)
        dI.sarj0 = vector(mode="numeric",length=n-1)
        dI.sar0k = vector(mode="numeric",length=n-1)
        dI.sbrj0 = vector(mode="numeric",length=n-1)
        dI.sbr0k = vector(mode="numeric",length=n-1)
        dI.sabjn = vector(mode="numeric",length=n-1)
        dI.sabnk = vector(mode="numeric",length=n-1)
        dI.sarjn = vector(mode="numeric",length=n-1)
        dI.sarnk = vector(mode="numeric",length=n-1)
        dI.sbrjn = vector(mode="numeric",length=n-1)
        dI.sbrnk = vector(mode="numeric",length=n-1)
        for (i in 1:(n-1)) {
            ## Midgut cells
            diffval.a <- get(paste0("I.ma",i-1)) - get(paste0("I.ma",i))
            diffval.b <- get(paste0("I.mb",i-1)) - get(paste0("I.mb",i))
            dI.ma[i] <- n * diffval.a / epsilon - beta.m * V.mb * get(paste0("I.ma",i))
            dI.mb[i] <- n * diffval.b / epsilon - beta.m * V.ma * get(paste0("I.mb",i))
            dI.mabj0[i] <- beta.m * V.mb * get(paste0("I.ma",i)) - n * get(paste0("I.mab",i,"0")) / epsilon
            dI.mab0k[i] <- beta.m * V.ma * get(paste0("I.mb",i)) - n * get(paste0("I.mab0",i)) / epsilon
            dI.mabjn[i] <- n * get(paste0("I.mab",i-1,n-1)) / epsilon
            dI.mabnk[i] <- n * get(paste0("I.mab",n-1,i-1)) / epsilon
            ## Secondary tissue virus
            sumV.ma.vec[i] <- get(paste0("I.mab",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.mabn",i)) * prod.fn.a(n,i,n)
            sumV.mb.vec[i] <- get(paste0("I.mab",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.mabn",i)) * prod.fn.b(n,i,n)
            sumV.mr.vec[i] <- (get(paste0("I.mab",i,"n")) + get(paste0("I.mabn",i))) * prod.fn.r(i,n,n)
            sumV.sa.vec[i] <- get(paste0("I.sab",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.sabn",i)) * prod.fn.a(n,i,n) +
                get(paste0("I.sar",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.sarn",i)) * prod.fn.a(n,i,n)
            sumV.sb.vec[i] <- get(paste0("I.sab",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.sabn",i)) * prod.fn.b(n,i,n) +
                get(paste0("I.sbr",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.sbrn",i)) * prod.fn.a(n,i,n)
            sumV.sr.vec[i] <- (get(paste0("I.sab",i,"n")) + get(paste0("I.sabn",i)) +
                               get(paste0("I.sar",i,"n")) + get(paste0("I.sarn",i)) +
                               get(paste0("I.sbr",i,"n")) + get(paste0("I.sbrn",i))) * prod.fn.r(i,n,n) +
                get(paste0("I.sar",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.sarn",i)) * prod.fn.b(n,i,n) +
                get(paste0("I.sbr",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.sbrn",i)) * prod.fn.b(n,i,n)
            ## Secondary tissue cells
            diffval.a <- get(paste0("I.sa",i-1)) - get(paste0("I.sa",i))
            diffval.b <- get(paste0("I.sb",i-1)) - get(paste0("I.sb",i))
            diffval.r <- get(paste0("I.sr",i-1)) - get(paste0("I.sr",i))
            dI.sa[i] <- n * diffval.a / tau - beta.s * (V.sb + V.sr) * get(paste0("I.sa",i))
            dI.sb[i] <- n * diffval.b / tau - beta.s * (V.sa + V.sr) * get(paste0("I.sb",i))
            dI.sr[i] <- n * diffval.r / tau - beta.s * (V.sa + V.sb) * get(paste0("I.sr",i))
            dI.sabj0[i] <- beta.s * V.sb * get(paste0("I.sa",i)) - n * get(paste0("I.sab",i,"0")) / tau
            dI.sab0k[i] <- beta.s * V.sa * get(paste0("I.sb",i)) - n * get(paste0("I.sab0",i)) / tau
            dI.sarj0[i] <- beta.s * V.sr * get(paste0("I.sa",i)) - n * get(paste0("I.sar",i,"0")) / tau
            dI.sar0k[i] <- beta.s * V.sa * get(paste0("I.sr",i)) - n * get(paste0("I.sar0",i)) / tau
            dI.sbrj0[i] <- beta.s * V.sr * get(paste0("I.sb",i)) - n * get(paste0("I.sbr",i,"0")) / tau
            dI.sbr0k[i] <- beta.s * V.sb * get(paste0("I.sr",i)) - n * get(paste0("I.sbr0",i)) / tau
            dI.sabjn[i] <- n * get(paste0("I.sab",i-1,n-1)) / tau
            dI.sabnk[i] <- n * get(paste0("I.sab",n-1,i-1)) / tau
            dI.sarjn[i] <- n * get(paste0("I.sar",i-1,n-1)) / tau
            dI.sarnk[i] <- n * get(paste0("I.sar",n-1,i-1)) / tau
            dI.sbrjn[i] <- n * get(paste0("I.sbr",i-1,n-1)) / tau
            dI.sbrnk[i] <- n * get(paste0("I.sbr",n-1,i-1)) / tau
        }
        ## Midgut cells
        dI.mabjk = array(dim = c(n-1,n-1))
        ## Secondary tissue cells
        dI.sabjk = array(dim = c(n-1,n-1))
        dI.sarjk = array(dim = c(n-1,n-1))
        dI.sbrjk = array(dim = c(n-1,n-1))
        for (j in 1:(n-1)) {
            for (k in 1:(n-1)) {
                ## Midgut cells
                diffval <- get(paste0("I.mab",j-1,k-1)) - get(paste0("I.mab",j,k))
                dI.mabjk[j,k] <- n * diffval / epsilon
                ## Secondary tissue cells
                diffval.ab <- get(paste0("I.sab",j-1,k-1)) - get(paste0("I.sab",j,k))
                dI.sabjk[j,k] <- n * diffval.ab / tau
                diffval.ar <- get(paste0("I.sar",j-1,k-1)) - get(paste0("I.sar",j,k))
                dI.sarjk[j,k] <- n * diffval.ar / tau
                diffval.br <- get(paste0("I.sbr",j-1,k-1)) - get(paste0("I.sbr",j,k))
                dI.sbrjk[j,k] <- n * diffval.br / tau                
            }
        }
        ## Midgut cells
        dI.mab00 <- beta.m * (V.ma * I.mb0 + V.mb * I.ma0) - n * I.mab00 / epsilon
        dI.mabnn <- n * get(paste0("I.mab",n-1,n-1)) / epsilon                
        ## Secondary tissue virus
        sumV.ma <- sum(sumV.ma.vec)
        sumV.mb <- sum(sumV.mb.vec)
        sumV.mr <- sum(sumV.mr.vec)
        sumV.sa <- sum(sumV.sa.vec)
        sumV.sb <- sum(sumV.sb.vec)
        sumV.sr <- sum(sumV.sr.vec)
        dV.sa = p.m * (I.man + sumV.ma + I.mabnn * prod.fn.a(n,n,n)) + p.s * (I.san + sumV.sa + (I.sabnn + I.sarnn) * prod.fn.a(n,n,n)) - c.s * V.sa
        dV.sb = p.m * (I.mbn + sumV.mb + I.mabnn * prod.fn.b(n,n,n)) + p.s * (I.sbn + sumV.sb + I.sabnn * prod.fn.b(n,n,n) + I.sbrnn * prod.fn.a(n,n,n)) - c.s * V.sb
        dV.sr = p.m * (sumV.mr + I.mabnn * prod.fn.r(n,n,n)) +
            p.s * (sumV.sr + (I.sabnn + I.sarnn + I.sbrnn) * prod.fn.r(n,n,n) + I.srn +
                   (I.sarnn + I.sbrnn) * prod.fn.b(n,n,n))  - c.s * V.sr
        ## Secondary tissue cells
        dI.sab00 <- beta.s * (V.sa * I.sb0 + V.sb * I.sa0) - n * I.sab00 / tau
        dI.sar00 <- beta.s * (V.sa * I.sr0 + V.sr * I.sa0) - n * I.sar00 / tau
        dI.sbr00 <- beta.s * (V.sr * I.sb0 + V.sb * I.sr0) - n * I.sbr00 / tau
        dI.sabnn <- n * get(paste0("I.sab",n-1,n-1)) / tau
        dI.sarnn <- n * get(paste0("I.sar",n-1,n-1)) / tau
        dI.sbrnn <- n * get(paste0("I.sbr",n-1,n-1)) / tau
        list(c(dV.ma, dV.mb,
               dT.m,
               dI.ma0, dI.ma, dI.man,
               dI.mb0, dI.mb, dI.mbn,
               dI.mab00, dI.mabj0, dI.mab0k,
               as.numeric(dI.mabjk),
               dI.mabnk, dI.mabjn, dI.mabnn,
               dV.sa, dV.sb, dV.sr,
               dT.s,
               dI.sa0, dI.sa, dI.san,
               dI.sb0, dI.sb, dI.sbn,
               dI.sr0, dI.sr, dI.srn,
               dI.sab00, dI.sabj0, dI.sab0k,
               as.numeric(dI.sabjk),
               dI.sabnk, dI.sabjn, dI.sabnn,
               dI.sar00, dI.sarj0, dI.sar0k,
               as.numeric(dI.sarjk),
               dI.sarnk, dI.sarjn, dI.sarnn,
               dI.sbr00, dI.sbrj0, dI.sbr0k,
               as.numeric(dI.sbrjk),
               dI.sbrnk, dI.sbrjn, dI.sbrnn))
    })
}

production.fn.a <- function(i,j,n,r=exp(n),wl){
    k <- r^(1/n)
    return(wl / (k^(j-i) + 1) + (1 - wl) / (k^(j-i) + 1)^10)
}
production.fn.b <- function(i,j,n,r=exp(n),wl){
    k <- r^(1/n)
    return(wl / (k^(i-j) + 1) + (1 - wl) / (k^(i-j) + 1)^10)
}
production.fn.r <- function(i,j,n,r=exp(n),wl){
    k <- r^(1/n)
    return((1 - wl)*(1 - (k^(10*(i-j)) + 1) / (k^(i-j) + 1)^10))
}




### Old functions, no longer use
## within midge model - original single infection model
wv.BTV.old = function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    dV.m = -c.m * V.m
    dT.m = -beta.m * T.m * V.m
    dI.m0 = beta.m * T.m * V.m - n * I.m0 / epsilon
    dI.m <- vector(mode="numeric",length=n-1)
    for (i in 1:(n-1)) {
        diffval <- get(paste0("I.m",i-1)) - get(paste0("I.m",i))
        dI.m[i] <- n * diffval / epsilon
    }
    dI.mn = n * get(paste0("I.m",n-1)) / epsilon
    ## dV.h = p * I.mn - c.s * V.h # only midgut infected cells produce virus
    dV.h = p.m * I.mn + p.s * I.sn - c.s * V.h # all infected cells produce virus
    dT.s = -beta.s * T.s * V.h
    dI.s0 = beta.s * T.s * V.h - n * I.s0 / tau
    dI.s <- vector(mode="numeric",length=n-1)
    for (i in 1:(n-1)) {
        diffval <- get(paste0("I.s",i-1)) - get(paste0("I.s",i))
        dI.s[i] <- n * diffval / tau
    }
    dI.sn = n * get(paste0("I.s",n-1)) / tau
    list(c(dV.m, dT.m, dI.m0, dI.m, dI.mn,
           dV.h, dT.s, dI.s0, dI.s, dI.sn))
  })
}

## within midge model intrathoracic - with eclipse
wv.BTV.intrathoracic.eclipse = function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    dV.h = p.s * I.sn - c.s * V.h # all infected cells produce virus
    dT.s = -beta.s * T.s * V.h
    dI.s0 = beta.s * T.s * V.h - n * I.s0 / tau
    dI.s <- vector(mode="numeric",length=n-1)
    if (n > 1){
        for (i in 1:(n-1)) {
            diffval <- get(paste0("I.s",i-1)) - get(paste0("I.s",i))
            dI.s[i] <- n * diffval / tau
        }
    }
    dI.sn = n * get(paste0("I.s",n-1)) / tau
    list(c(dV.h, dT.s, dI.s0, dI.s, dI.sn))
  })
}

## within midge coinfection
wv.BTV.coinfection = function(t, state, parameters, prod.fn.a=production.fn.a,
                              prod.fn.b=production.fn.b, prod.fn.r=production.fn.r)
{
    with(as.list(c(state, parameters)),{
        ## Midgut virus
        dV.ma = -c.m * V.ma
        dV.mb = -c.m * V.mb
        ## Midgut cells
        dT.m = -beta.m * T.m * (V.ma + V.mb)
        dI.ma0 = beta.m * T.m * V.ma - I.ma0 * (n / epsilon + beta * V.mb)
        dI.mb0 = beta.m * T.m * V.mb - I.mb0 * (n / epsilon + beta * V.ma)
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
        sumV.sa.vec <- vector(mode="numeric",length=n-1)
        sumV.sb.vec <- vector(mode="numeric",length=n-1)
        sumV.sr.vec <- vector(mode="numeric",length=n-1)
        ## Secondary tissue cells
        dT.s = -beta.s * T.s * (V.sa + V.sb)
        dI.sa0 = beta.s * T.s * V.sa - I.sa0 * (n / tau + beta * V.sb)
        dI.sb0 = beta.s * T.s * V.sb - I.sb0 * (n / tau + beta * V.sa)
        dI.san = n * get(paste0("I.sa",n-1)) / tau
        dI.sbn = n * get(paste0("I.sb",n-1)) / tau
        dI.sa <- vector(mode="numeric",length=n-1)
        dI.sb <- vector(mode="numeric",length=n-1)
        dI.sabj0 = vector(mode="numeric",length=n-1)
        dI.sab0k = vector(mode="numeric",length=n-1)
        dI.sabjn = vector(mode="numeric",length=n-1)
        dI.sabnk = vector(mode="numeric",length=n-1)
        for (i in 1:(n-1)) {
            ## Midgut cells
            diffval.a <- get(paste0("I.ma",i-1)) - get(paste0("I.ma",i))
            diffval.b <- get(paste0("I.mb",i-1)) - get(paste0("I.mb",i))
            dI.ma[i] <- n * diffval.a / epsilon - beta.m * V.mb * get(paste0("I.ma",i))
            dI.mb[i] <- n * diffval.b / epsilon - beta.m * V.ma * get(paste0("I.mb",i))
            dI.mabj0[i] <- beta.m * V.mb * get(paste0("I.ma",i)) - n * get(paste0("I.mab",i,"0")) / epsilon
            dI.mab0k[i] <- beta.m * V.ma * get(paste0("I.mb",i)) - n * get(paste0("I.mab0",i)) / epsilon
            dI.mabjn[i] <- n * get(paste0("I.mab",i-1,n-1)) / epsilon
            dI.mabnk[i] <- n * get(paste0("I.mab",n-1,i-1)) / epsilon
            ## Secondary tissue virus
            sumV.ma.vec[i] <- get(paste0("I.mab",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.mabn",i)) * prod.fn.a(n,i,n)
            sumV.mb.vec[i] <- get(paste0("I.mab",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.mabn",i)) * prod.fn.b(n,i,n)
            sumV.mr.vec[i] <- (get(paste0("I.mab",i,"n")) + get(paste0("I.mabn",i))) * prod.fn.r(i,n,n)
            sumV.sa.vec[i] <- get(paste0("I.sab",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.sabn",i)) * prod.fn.a(n,i,n)
            sumV.sb.vec[i] <- get(paste0("I.sab",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.sabn",i)) * prod.fn.b(n,i,n)
            sumV.sr.vec[i] <- (get(paste0("I.sab",i,"n")) + get(paste0("I.sabn",i))) * prod.fn.r(i,n,n)
            ## Secondary tissue cells
            diffval.a <- get(paste0("I.sa",i-1)) - get(paste0("I.sa",i))
            diffval.b <- get(paste0("I.sb",i-1)) - get(paste0("I.sb",i))
            dI.sa[i] <- n * diffval.a / tau - beta.s * V.sb * get(paste0("I.sa",i))
            dI.sb[i] <- n * diffval.b / tau - beta.s * V.sa * get(paste0("I.sb",i))
            dI.sabj0[i] <- beta.s * V.sb * get(paste0("I.sa",i)) - n * get(paste0("I.sab",i,"0")) / tau
            dI.sab0k[i] <- beta.s * V.sa * get(paste0("I.sb",i)) - n * get(paste0("I.sab0",i)) / tau
            dI.sabjn[i] <- n * get(paste0("I.sab",i-1,n-1)) / tau
            dI.sabnk[i] <- n * get(paste0("I.sab",n-1,i-1)) / tau
        }
        ## Midgut cells
        dI.mabjk = array(dim = c(n-1,n-1))
        ## Secondary tissue cells
        dI.sabjk = array(dim = c(n-1,n-1))
        for (j in 1:(n-1)) {
            for (k in 1:(n-1)) {
                ## Midgut cells
                diffval <- get(paste0("I.mab",j-1,k-1)) - get(paste0("I.mab",j,k))
                dI.mabjk[j,k] <- n * diffval / epsilon
                ## Secondary tissue cells
                diffval <- get(paste0("I.sab",j-1,k-1)) - get(paste0("I.sab",j,k))
                dI.sabjk[j,k] <- n * diffval / tau                
            }
        }
        ## Midgut cells
        dI.mab00 <- beta.m * (V.ma * I.mb0 + V.mb * I.ma0) - n * I.mab00 / epsilon
        dI.mcn <- n * get(paste0("I.mab",n-1,n-1)) / epsilon                
        ## Secondary tissue virus
        sumV.ma <- sum(sumV.ma.vec)
        sumV.mb <- sum(sumV.mb.vec)
        sumV.mr <- sum(sumV.mr.vec)
        sumV.sa <- sum(sumV.sa.vec)
        sumV.sb <- sum(sumV.sb.vec)
        sumV.sr <- sum(sumV.sr.vec)
        dV.sa = p.m * (I.man + sumV.ma + I.mcn * prod.fn.a(n,n,n)) + p.s * (I.san + sumV.sa + I.scn * prod.fn.a(n,n,n)) - c.s * V.sa
        dV.sb = p.m * (I.mbn + sumV.mb + I.mcn * prod.fn.b(n,n,n)) + p.s * (I.sbn + sumV.sb + I.scn * prod.fn.b(n,n,n)) - c.s * V.sb
        dV.sr = p.m * (sumV.mr + I.mcn * prod.fn.r(n,n,n)) + p.s * (sumV.sr + I.scn * prod.fn.r(n,n,n)) - c.s * V.sr
        ## Secondary tissue cells
        dI.sab00 <- beta.s * (V.sa * I.sb0 + V.sb * I.sa0) - n * I.sab00 / tau
        dI.scn <- n * get(paste0("I.sab",n-1,n-1)) / tau
        list(c(dV.ma, dV.mb,
               dT.m,
               dI.ma0, dI.ma, dI.man,
               dI.mb0, dI.mb, dI.mbn,
               dI.mab00, dI.mabj0, dI.mab0k,
               as.numeric(dI.mabjk),
               dI.mabnk, dI.mabjn, dI.mcn,
               dV.sa, dV.sb, dV.sr,
               dT.s,
               dI.sa0, dI.sa, dI.san,
               dI.sb0, dI.sb, dI.sbn,
               dI.sab00, dI.sabj0, dI.sab0k,
               as.numeric(dI.sabjk),
               dI.sabnk, dI.sabjn, dI.scn))
    })
}

## within midge reassortment NEEDS EDITING
wv.BTV.coinfection.reassort.old = function(t, state, parameters, prod.fn.a=production.fn.a,
                                       prod.fn.b=production.fn.b, prod.fn.r=production.fn.r)
{
    with(as.list(c(state, parameters)),{
        ## Midgut virus
        dV.ma = -c.m * V.ma
        dV.mb = -c.m * V.mb
        ## Midgut cells
        dT.m = -beta.m * T.m * (V.ma + V.mb)
        dI.ma0 = beta.m * T.m * V.ma - I.ma0 * (n / epsilon + beta * V.mb)
        dI.mb0 = beta.m * T.m * V.mb - I.mb0 * (n / epsilon + beta * V.ma)
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
        sumV.sa.vec <- vector(mode="numeric",length=n-1)
        sumV.sb.vec <- vector(mode="numeric",length=n-1)
        sumV.sr.vec <- vector(mode="numeric",length=n-1)
        ## Secondary tissue cells
        dT.s = -beta.s * T.s * (V.sa + V.sb + V.sr)
        dI.sa0 = beta.s * T.s * V.sa - I.sa0 * (n / tau + beta * (V.sb + V.sr))
        dI.sb0 = beta.s * T.s * V.sb - I.sb0 * (n / tau + beta * (V.sa + V.sr))
        dI.sr0 = beta.s * T.s * V.sr - I.sr0 * (n / tau + beta * (V.sa + V.sb))
        dI.san = n * get(paste0("I.sa",n-1)) / tau
        dI.sbn = n * get(paste0("I.sb",n-1)) / tau
        dI.srn = n * get(paste0("I.sr",n-1)) / tau
        dI.sa <- vector(mode="numeric",length=n-1)
        dI.sb <- vector(mode="numeric",length=n-1)
        dI.sr <- vector(mode="numeric",length=n-1)
        dI.sabj0 = vector(mode="numeric",length=n-1)
        dI.sab0k = vector(mode="numeric",length=n-1)
        dI.sarj0 = vector(mode="numeric",length=n-1)
        dI.sar0k = vector(mode="numeric",length=n-1)
        dI.sbrj0 = vector(mode="numeric",length=n-1)
        dI.sbr0k = vector(mode="numeric",length=n-1)
        dI.sabjn = vector(mode="numeric",length=n-1)
        dI.sabnk = vector(mode="numeric",length=n-1)
        dI.sarjn = vector(mode="numeric",length=n-1)
        dI.sarnk = vector(mode="numeric",length=n-1)
        dI.sbrjn = vector(mode="numeric",length=n-1)
        dI.sbrnk = vector(mode="numeric",length=n-1)
        for (i in 1:(n-1)) {
            ## Midgut cells
            diffval.a <- get(paste0("I.ma",i-1)) - get(paste0("I.ma",i))
            diffval.b <- get(paste0("I.mb",i-1)) - get(paste0("I.mb",i))
            dI.ma[i] <- n * diffval.a / epsilon - beta.m * V.mb * get(paste0("I.ma",i))
            dI.mb[i] <- n * diffval.b / epsilon - beta.m * V.ma * get(paste0("I.mb",i))
            dI.mabj0[i] <- beta.m * V.mb * get(paste0("I.ma",i)) - n * get(paste0("I.mab",i,"0")) / epsilon
            dI.mab0k[i] <- beta.m * V.ma * get(paste0("I.mb",i)) - n * get(paste0("I.mab0",i)) / epsilon
            dI.mabjn[i] <- n * get(paste0("I.mab",i-1,n-1)) / epsilon
            dI.mabnk[i] <- n * get(paste0("I.mab",n-1,i-1)) / epsilon
            ## Secondary tissue virus
            sumV.ma.vec[i] <- get(paste0("I.mab",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.mabn",i)) * prod.fn.a(n,i,n)
            sumV.mb.vec[i] <- get(paste0("I.mab",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.mabn",i)) * prod.fn.b(n,i,n)
            sumV.mr.vec[i] <- (get(paste0("I.mab",i,"n")) + get(paste0("I.mabn",i))) * prod.fn.r(i,n,n)
            sumV.sa.vec[i] <- get(paste0("I.sab",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.sabn",i)) * prod.fn.a(n,i,n) +
                get(paste0("I.sar",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.sarn",i)) * prod.fn.a(n,i,n)
            sumV.sb.vec[i] <- get(paste0("I.sab",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.sabn",i)) * prod.fn.b(n,i,n) +
                get(paste0("I.sbr",i,"n")) * prod.fn.a(i,n,n) + get(paste0("I.sbrn",i)) * prod.fn.a(n,i,n)
            sumV.sr.vec[i] <- (get(paste0("I.sab",i,"n")) + get(paste0("I.sabn",i)) +
                               get(paste0("I.sar",i,"n")) + get(paste0("I.sarn",i)) +
                               get(paste0("I.sbr",i,"n")) + get(paste0("I.sbrn",i))) * prod.fn.r(i,n,n) +
                get(paste0("I.sar",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.sarn",i)) * prod.fn.b(n,i,n) +
                get(paste0("I.sbr",i,"n")) * prod.fn.b(i,n,n) + get(paste0("I.sbrn",i)) * prod.fn.b(n,i,n)
            ## Secondary tissue cells
            diffval.a <- get(paste0("I.sa",i-1)) - get(paste0("I.sa",i))
            diffval.b <- get(paste0("I.sb",i-1)) - get(paste0("I.sb",i))
            diffval.r <- get(paste0("I.sr",i-1)) - get(paste0("I.sr",i))
            dI.sa[i] <- n * diffval.a / tau - beta.s * (V.sb + V.sr) * get(paste0("I.sa",i))
            dI.sb[i] <- n * diffval.b / tau - beta.s * (V.sa + V.sr) * get(paste0("I.sb",i))
            dI.sr[i] <- n * diffval.r / tau - beta.s * (V.sa + V.sb) * get(paste0("I.sr",i))
            dI.sabj0[i] <- beta.s * V.sb * get(paste0("I.sa",i)) - n * get(paste0("I.sab",i,"0")) / tau
            dI.sab0k[i] <- beta.s * V.sa * get(paste0("I.sb",i)) - n * get(paste0("I.sab0",i)) / tau
            dI.sarj0[i] <- beta.s * V.sr * get(paste0("I.sa",i)) - n * get(paste0("I.sar",i,"0")) / tau
            dI.sar0k[i] <- beta.s * V.sa * get(paste0("I.sr",i)) - n * get(paste0("I.sar0",i)) / tau
            dI.sbrj0[i] <- beta.s * V.sr * get(paste0("I.sb",i)) - n * get(paste0("I.sbr",i,"0")) / tau
            dI.sbr0k[i] <- beta.s * V.sb * get(paste0("I.sr",i)) - n * get(paste0("I.sbr0",i)) / tau
            dI.sabjn[i] <- n * get(paste0("I.sab",i-1,n-1)) / tau
            dI.sabnk[i] <- n * get(paste0("I.sab",n-1,i-1)) / tau
            dI.sarjn[i] <- n * get(paste0("I.sar",i-1,n-1)) / tau
            dI.sarnk[i] <- n * get(paste0("I.sar",n-1,i-1)) / tau
            dI.sbrjn[i] <- n * get(paste0("I.sbr",i-1,n-1)) / tau
            dI.sbrnk[i] <- n * get(paste0("I.sbr",n-1,i-1)) / tau
        }
        ## Midgut cells
        dI.mabjk = array(dim = c(n-1,n-1))
        ## Secondary tissue cells
        dI.sabjk = array(dim = c(n-1,n-1))
        dI.sarjk = array(dim = c(n-1,n-1))
        dI.sbrjk = array(dim = c(n-1,n-1))
        for (j in 1:(n-1)) {
            for (k in 1:(n-1)) {
                ## Midgut cells
                diffval <- get(paste0("I.mab",j-1,k-1)) - get(paste0("I.mab",j,k))
                dI.mabjk[j,k] <- n * diffval / epsilon
                ## Secondary tissue cells
                diffval.ab <- get(paste0("I.sab",j-1,k-1)) - get(paste0("I.sab",j,k))
                dI.sabjk[j,k] <- n * diffval.ab / tau
                diffval.ar <- get(paste0("I.sar",j-1,k-1)) - get(paste0("I.sar",j,k))
                dI.sarjk[j,k] <- n * diffval.ar / tau
                diffval.br <- get(paste0("I.sbr",j-1,k-1)) - get(paste0("I.sbr",j,k))
                dI.sbrjk[j,k] <- n * diffval.br / tau                
            }
        }
        ## Midgut cells
        dI.mab00 <- beta.m * (V.ma * I.mb0 + V.mb * I.ma0) - n * I.mab00 / epsilon
        dI.mabnn <- n * get(paste0("I.mab",n-1,n-1)) / epsilon                
        ## Secondary tissue virus
        sumV.ma <- sum(sumV.ma.vec)
        sumV.mb <- sum(sumV.mb.vec)
        sumV.mr <- sum(sumV.mr.vec)
        sumV.sa <- sum(sumV.sa.vec)
        sumV.sb <- sum(sumV.sb.vec)
        sumV.sr <- sum(sumV.sr.vec)
        dV.sa = p.m * (I.man + sumV.ma + I.mabnn * prod.fn.a(n,n,n)) + p.s * (I.san + sumV.sa + (I.sabnn + I.sarnn) * prod.fn.a(n,n,n)) - c.s * V.sa
        dV.sb = p.m * (I.mbn + sumV.mb + I.mabnn * prod.fn.b(n,n,n)) + p.s * (I.sbn + sumV.sb + I.sabnn * prod.fn.b(n,n,n) + I.sbrnn * prod.fn.a(n,n,n)) - c.s * V.sb
        dV.sr = p.m * (sumV.mr + I.mabnn * prod.fn.r(n,n,n)) +
            p.s * (sumV.sr + (I.sabnn + I.sarnn + I.sbrnn) * prod.fn.r(n,n,n) + I.srn +
                   (I.sarnn + I.sbrnn) * prod.fn.b(n,n,n))  - c.s * V.sr
        ## Secondary tissue cells
        dI.sab00 <- beta.s * (V.sa * I.sb0 + V.sb * I.sa0) - n * I.sab00 / tau
        dI.sar00 <- beta.s * (V.sa * I.sr0 + V.sr * I.sa0) - n * I.sar00 / tau
        dI.sbr00 <- beta.s * (V.sr * I.sb0 + V.sb * I.sr0) - n * I.sbr00 / tau
        dI.sabnn <- n * get(paste0("I.sab",n-1,n-1)) / tau
        dI.sarnn <- n * get(paste0("I.sar",n-1,n-1)) / tau
        dI.sbrnn <- n * get(paste0("I.sbr",n-1,n-1)) / tau
        list(c(dV.ma, dV.mb,
               dT.m,
               dI.ma0, dI.ma, dI.man,
               dI.mb0, dI.mb, dI.mbn,
               dI.mab00, dI.mabj0, dI.mab0k,
               as.numeric(dI.mabjk),
               dI.mabnk, dI.mabjn, dI.mabnn,
               dV.sa, dV.sb, dV.sr,
               dT.s,
               dI.sa0, dI.sa, dI.san,
               dI.sb0, dI.sb, dI.sbn,
               dI.sr0, dI.sr, dI.srn,
               dI.sab00, dI.sabj0, dI.sab0k,
               as.numeric(dI.sabjk),
               dI.sabnk, dI.sabjn, dI.sabnn,
               dI.sar00, dI.sarj0, dI.sar0k,
               as.numeric(dI.sarjk),
               dI.sarnk, dI.sarjn, dI.sarnn,
               dI.sbr00, dI.sbrj0, dI.sbr0k,
               as.numeric(dI.sbrjk),
               dI.sbrnk, dI.sbrjn, dI.sbrnn))
    })
}
