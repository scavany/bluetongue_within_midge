# set working directory
setwd('~/Documents/Research /')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(deSolve)){install.packages('deSolve'); library(deSolve)}

# specify within-vector model of BTV infection
wv.BTV = function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    dV.m = -beta.me * T.me * V.m - k.m * V.m
    dT.me = -beta.me * T.me * V.m
    dE.me.1 = beta.me * T.me * V.m - alpha * n * E.me.1
    dE.me.2 = alpha * n * (E.me.1 - E.me.2)
    dI.me.2 = alpha * n * E.me.2 - k.me * I.me.2
    dV.h = p.me * I.me.2 - k.h * V.h - beta.st * T.st * V.h
    dT.st = -beta.st * T.st * V.h
    dE.st.1 = beta.st * T.st * V.h - alpha * n * E.st.1
    dE.st.2  = alpha * n * (E.st.1 - E.st.2)
    dI.st.2 = alpha * n * E.st.2 - k.st * I.st.2
    dV.s = p.st * I.st.2 - k.s * V.s 
    list(c(dV.m, dT.me, dE.me.1, dE.me.2, dI.me.2, dV.h,
           dT.st, dE.st.1, dE.st.2, dI.st.2, dV.s))
  })
}


# run model
times = seq(from = 0, to = 50, by = 0.01)
parameters = c(beta.me = 0.001, k.m = 2, n = 2, alpha = 10, p.me = 10, k.h = 1,
               beta.st = 0.1, epsilon = 0.2, p.st = 10, k.s = 1, k.me = 0.2, k.st = 0.2)
state = c(V.m = 100, T.me = 100, E.me.1 = 0, E.me.2 = 0, I.me.2 = 0, V.h = 0,
          T.st = 100, E.st.1 = 0, E.st.2 = 0, I.st.2 = 0, V.s = 0)

dat = ode(y = state, times = times, func = wv.BTV, parms = parameters)
1 * ((dat[,2] + dat[,7] + dat[,12]) > 300)

# generate plot
par(mar=c(5.1,4.1,4.1,2.1))
layout(mat = matrix(1:8, nrow = 2))
plot(dat[,1], dat[,2], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Midgut Viral Load', las = 1)
plot(dat[,1], dat[,3], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'Uninfected Midgut Epithelial Cells', las = 1)
plot(dat[,1], dat[,4] + dat[,5] + dat[,6], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Infected Midgut Epithelial Cells', las = 1)
plot(dat[,1], dat[,7], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'Hemocoel Viral Load', las = 1)
plot(dat[,1], dat[,8], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'Uninfected Secondary Tissue Cells', las = 1)
plot(dat[,1], dat[,9] + dat[,10] + dat[,11], type = 'l', lwd = 2, bty = 'n',
     xlab = 'Time (Days)', ylab = 'Infected Secondary Tissue Cells', las = 1)
plot(dat[,1], dat[,12], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Salivary Gland Viral Load', las = 1)
plot(dat[,1], dat[,2] + dat[,7] + dat[,12], type = 'l', lwd = 2, bty = 'n', 
     xlab = 'Time (Days)', ylab = 'Total Viral Load', las = 1, ylim = c(0,800))

