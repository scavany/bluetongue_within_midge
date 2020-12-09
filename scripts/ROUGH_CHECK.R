## Secondary tissue virus
sumV.ma.vec <- vector(mode="numeric",length=n-1)
sumV.mb.vec <- vector(mode="numeric",length=n-1)
sumV.mr.vec <- vector(mode="numeric",length=n-1)
sumV.sa.vec <- vector(mode="numeric",length=n-1)
sumV.sb.vec <- vector(mode="numeric",length=n-1)
sumV.sr.vec <- vector(mode="numeric",length=n-1)
## Secondary tissue cells
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
    dI.sabj0[i] <- - n * get(paste0("I.sab",i,"0")) / tau
    dI.sab0k[i] <- - n * get(paste0("I.sab0",i)) / tau
    dI.sarj0[i] <- - n * get(paste0("I.sar",i,"0")) / tau
    dI.sar0k[i] <- - n * get(paste0("I.sar0",i)) / tau
    dI.sbrj0[i] <- - n * get(paste0("I.sbr",i,"0")) / tau
    dI.sbr0k[i] <- - n * get(paste0("I.sbr0",i)) / tau
    dI.sabjn[i] <- n * get(paste0("I.sab",i-1,n-1)) / tau
    dI.sabnk[i] <- n * get(paste0("I.sab",n-1,i-1)) / tau
    dI.sarjn[i] <- n * get(paste0("I.sar",i-1,n-1)) / tau
    dI.sarnk[i] <- n * get(paste0("I.sar",n-1,i-1)) / tau
    dI.sbrjn[i] <- n * get(paste0("I.sbr",i-1,n-1)) / tau
    dI.sbrnk[i] <- n * get(paste0("I.sbr",n-1,i-1)) / tau
}
## Secondary tissue cells
dI.sabjk = array(dim = c(n-1,n-1))
dI.sarjk = array(dim = c(n-1,n-1))
dI.sbrjk = array(dim = c(n-1,n-1))
for (j in 1:(n-1)) {
    for (k in 1:(n-1)) {
        ## Secondary tissue cells
        diffval.ab <- get(paste0("I.sab",j-1,k-1)) - get(paste0("I.sab",j,k))
        dI.sabjk[j,k] <- n * diffval.ab / tau
        diffval.ar <- get(paste0("I.sar",j-1,k-1)) - get(paste0("I.sar",j,k))
        dI.sarjk[j,k] <- n * diffval.ar / tau
        diffval.br <- get(paste0("I.sbr",j-1,k-1)) - get(paste0("I.sbr",j,k))
        dI.sbrjk[j,k] <- n * diffval.br / tau                
    }
}
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
dI.sab00 <- - n * I.sab00 / tau
dI.sar00 <- - n * I.sar00 / tau
dI.sbr00 <- - n * I.sbr00 / tau
dI.sabnn <- n * get(paste0("I.sab",n-1,n-1)) / tau
dI.sarnn <- n * get(paste0("I.sar",n-1,n-1)) / tau
dI.sbrnn <- n * get(paste0("I.sbr",n-1,n-1)) / tau
