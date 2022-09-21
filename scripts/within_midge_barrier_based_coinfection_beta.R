## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspacer
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(data.table,
       viridis,
       deSolve,
       grDevices,
       BayesianTools,
       mixtools,
       adaptivetau,
       mgcv,
       RColorBrewer,
       imager,
       sensobol#,psych,bbmle
       )

## Plot controls
check.neutrality <- FALSE # check for eco. and pop gen neutrality (Lipsitch et al)
check.works <- FALSE  # check coinfection behaves the same way as the single model
plot.timeseries <- FALSE # plot coinfection timeseries
plot.reassortment <- FALSE # plot reassortment proportion by gap
regenerate.parmsweep <- FALSE
plot.parmsweep <- TRUE
plot.sobol <- FALSE

## load within-vector model of BTV infection
source("./within_midge_barrier_based_fn.R")

## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
initial.titre <- 10^fu$log10.mean.titre[1]
equilibrium.titre <- mean(tail(10^fu$log10.mean.titre)) ##Should this be the mean of the log10??
lod <- 10^0.75

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)
initial.titre.it <- fu.intrathoracic$titre[1]
equilibrium.titre.it <- mean(tail(fu.intrathoracic$titre)) ##Should this be the mean of the log10??

## Set up timing vectors
times = seq(from = 0, to = 384, by = 1)
data.hours.it <- which(times %in% fu.intrathoracic.hours)
model.hours.it <- which(fu.intrathoracic.hours %in% times)
titre.data.it <- floor(fu.intrathoracic$titre[model.hours.it]+0.5)
data.hours <- which(times %in% fu$Time.pi..hours.)
model.hours <- which(fu$Time.pi..hours. %in% times)
titre.data <- floor(10^fu$log10.mean.titre[model.hours] + 0.5)

## Fit logistic curve to the proportion positive
prop.positive.fn <- function(times,prop.disseminated,k,t0) {
    1 - (1 - prop.disseminated)/(1 + exp(-k*(times - t0)))
}
logistic.optim.fn <- function(par) {
    prop.positive <- prop.positive.fn(times,par[1],par[2],par[3])
    return(sum((prop.positive[data.hours] - fu$Detection.rate[model.hours])^2))
}
logistic.optim.out <- optim(c(0.4,1,0),logistic.optim.fn)
prop.positive.vec <- prop.positive.fn(times,logistic.optim.out$par[1],
                                      logistic.optim.out$par[2],
                                      logistic.optim.out$par[3])
prop.disseminated <- min(prop.positive.vec)

## Get expression for cb and cm
c.m <- log(initial.titre/lod) / logistic.optim.out$par[3]

## Barrier probabilities
p.mib=1/3;p.db=0.12
barrier.probs <- c(p.mib=p.mib,p.db=p.db)

### Load mcmc results
##load(file="./barrier_mcmc_out.RData",verbose=T) # fit to all data
##load(file="./barrier_mcmc_out_oral_only.RData",verbose=T) # fit only to oral data
##load(file="./barrier_mcmc_out_backflow.RData",verbose=T) # allow backflow
##load(file="./barrier_mcmc_out_backflow_oral_only.RData",verbose=T) # fit to oral only, allow backflow
load(file="./coinfection_mcmc_out_joint_fitting_10000.RData",verbose=T) # fit to all data

sample.out <- as.data.frame(getSample(mcmc.out,start=mcmc.out$settings$iterations / 9))
MAP.out <- MAP(mcmc.out,start=mcmc.out$settings$iterations / 9)$parametersMAP
CI.out <- apply(sample.out,2,function(x)quantile(x,c(0.025,0.975)))

## Now get new parameters and plot
parms <- MAP.out
parms["beta.m"] <- parms["beta.mT.m"]/parms["T0.m"]
parms["c.m"] <- c.m;parms["c.b"] <- c.m;parms["beta.b"] <- parms["beta.m"];parms["mui.m"] <- parms["mu.m"];parms["mui.s"] <- parms["mu.s"]
state <- c(V.b = floor(initial.titre*log(2) + 0.5),
           V.m=0, T.m=parms[["T0.m"]], E.m = 0, E.mm = 0, I.m=0, I.mm=0, 
           V.s=0, T.s=parms[["T0.s"]], E.s = 0, E.ss = 0, I.s=0, I.ss=0
           )
state.coinf <- c(V.ba = floor(initial.titre*log(2) + 0.5), V.bb = 0, V.br = 0,
                 V.ma = 0, V.mb = 0, V.mr = 0,
                 T.m = parms[["T0.m"]],
                 E.ma = 0,E.mb = 0,E.mr = 0,
                 E.maa = 0,E.mbb = 0,E.mrr = 0,E.mab = 0,E.mar = 0,E.mbr = 0,
                 I.ma = 0,I.mb = 0,I.mr = 0,
                 I.maa = 0,I.mbb = 0,I.mrr = 0,I.mab = 0,I.mar = 0,I.mbr = 0,
                 V.sa = 0, V.sb = 0, V.sr = 0,
                 T.s = parms[["T0.s"]],
                 E.sa = 0,E.sb = 0,E.sr = 0,
                 E.saa = 0,E.sbb = 0,E.srr = 0,E.sab = 0,E.sar = 0,E.sbr = 0,
                 I.sa = 0,I.sb = 0,I.sr = 0,
                 I.saa = 0,I.sbb = 0,I.srr = 0,I.sab = 0,I.sar = 0,I.sbr = 0)

## Check for population genetic neutrality in reduced model
if (check.neutrality) {
    state.coinf.simple <- c(V.ba = floor(initial.titre*log(2) + 0.5), V.bb = 0,
                            V.ma = 0, V.mb = 0,
                            T.m = parms[["T0.m"]],
                            E.ma = 0,E.mb = 0,E.maa = 0,E.mbb = 0,E.mab = 0,
                            I.ma = 0,I.mb = 0,I.maa = 0,I.mbb = 0,I.mab = 0)
    state.coinf.simple.50 <- state.coinf.simple
    state.coinf.simple.50["V.ba"] <- 0.5 * floor(initial.titre*log(2) + 0.5)
    state.coinf.simple.50["V.bb"] <- 0.5 * floor(initial.titre*log(2) + 0.5)
    coinfection.out.50 <- wv.BTV.coinfection.simple(times,state.coinf.simple.50,parms)
    totalA.50 <- coinfection.out.50$V.ma + coinfection.out.50$V.ba
    totalB.50 <- coinfection.out.50$V.mb + coinfection.out.50$V.bb
    state.coinf.simple.75 <- state.coinf.simple
    state.coinf.simple.75["V.ba"] <- 0.75 * floor(initial.titre*log(2) + 0.5)
    state.coinf.simple.75["V.bb"] <- 0.25 * floor(initial.titre*log(2) + 0.5)
    coinfection.out.75 <- wv.BTV.coinfection.simple(times,state.coinf.simple.75,parms)
    totalA.75 <- coinfection.out.75$V.ma + coinfection.out.75$V.ba
    totalB.75 <- coinfection.out.75$V.mb + coinfection.out.75$V.bb
    state.coinf.simple.90 <- state.coinf.simple
    state.coinf.simple.90["V.ba"] <- 0.90 * floor(initial.titre*log(2) + 0.5)
    state.coinf.simple.90["V.bb"] <- 0.10 * floor(initial.titre*log(2) + 0.5)
    coinfection.out.90 <- wv.BTV.coinfection.simple(times,state.coinf.simple.90,parms)
    totalA.90 <- coinfection.out.90$V.ma + coinfection.out.90$V.ba
    totalB.90 <- coinfection.out.90$V.mb + coinfection.out.90$V.bb
    
    cols <- brewer.pal(4,"Reds")[2:4]
    cols.B <- brewer.pal(4,"Blues")[2:4]
    pdf("../figures/popgen_neutrality.pdf")
    plot(coinfection.out.50$time,
         totalA.50/rowSums(coinfection.out.50[,grepl("V.",names(coinfection.out.50))]),
         ylim = c(0,1), col=cols[1],lwd=3,type='l',xlab="Time (hours)",ylab="Genotype frequency")
    lines(coinfection.out.50$time,
          totalB.50/rowSums(coinfection.out.50[,grepl("V.",names(coinfection.out.90))]),
          col=cols.B[1],lwd=3,lty=2)     
    lines(coinfection.out.75$time,
          totalA.75/rowSums(coinfection.out.75[,grepl("V.",names(coinfection.out.75))]),
          col=cols[2],lwd=3)
    lines(coinfection.out.75$time,
          totalB.75/rowSums(coinfection.out.75[,grepl("V.",names(coinfection.out.75))]),
          col=cols.B[2],lwd=3,lty=2)
    lines(coinfection.out.90$time,
          totalA.90/rowSums(coinfection.out.90[,grepl("V.",names(coinfection.out.90))]),
          col=cols[3],lwd=3)     
    lines(coinfection.out.90$time,
          totalB.90/rowSums(coinfection.out.90[,grepl("V.",names(coinfection.out.90))]),
          col=cols.B[3],lwd=3,lty=2)
    text(330,0.68,"Initial frequency of A")
    legend(267,0.68,legend=c("","",""),
           col=cols,lty=1,lwd=3,bty="n")
    legend(297,0.68,legend=c("50%","75%","90%"),
           col=cols.B,lty=1,lwd=3,bty="n")
    legend(267,0.44,legend=c("Genotype A","Genotype B"),
           lty=c(1,2),lwd=3,bty="n")
    dev.off()
}

## Check coinfection model works similarly for single infection
wl <- 1
cols <- rev(viridis(4))
if (check.works) {
    pdf("../figures/model_equivalence.pdf")
    single.out <- wv.BTV.barrier.det(times,state,parms)
    V.out.single <- rowSums(single.out[,grepl("V.",names(single.out))])
    coinfection.out <- wv.BTV.coinfection.reassort(times,state.coinf,parms,withlike=wl)
    V.out.a <- rowSums(coinfection.out[,grepl("V.",names(coinfection.out))])
    plot(times,V.out.a,
         type='l',log='y',lwd=9,col=cols[1],
         las=1,xlab="Time since infection (hours)",
         ylab=expression("Viral titre ("~TCID[50]~")"),
         bty="n")
    lines(times,rowSums(single.out[,grepl("V.",names(single.out))]),
          lty=1,col=cols[2],lwd=3)
    print(paste("A vs single",
                sum(abs(V.out.a - V.out.single))))
    state.coinf[["V.ba"]] <- 0; state.coinf[["V.bb"]] <- floor(initial.titre*log(2) + 0.5);
    coinfection.out <- wv.BTV.coinfection.reassort(times,state.coinf,parms,withlike=wl)
    V.out.b <- rowSums(coinfection.out[,grepl("V.",names(coinfection.out))])
    print(paste("A vs B",
                sum(abs(V.out.a - V.out.b))))
    lines(times,V.out.b,
          col=cols[3],lwd=3,lty=2)
    state.coinf[["V.ba"]] <- 0.5 * floor(initial.titre*log(2) + 0.5)
    state.coinf[["V.bb"]] <- 0.5 * floor(initial.titre*log(2) + 0.5)
    coinfection.out <- wv.BTV.coinfection.reassort(times,state.coinf,parms,withlike=wl)
    V.out.mix <- rowSums(coinfection.out[,grepl("V.",names(coinfection.out))])
    lines(times,V.out.mix,
          col=cols[4],lwd=3,lty=3)
    print(paste("A vs mixture",
                sum(abs(V.out.a - V.out.mix))))
    legend("right",bty="n",legend=c("Single infection model","Coinfection model, just A",
                                    "Coinfection model, just B", "Coinfection model, equal mix"),
           col=cols,lty=c(1,1:3),lwd=c(9,rep(3,3)))
    dev.off()
}

## Samal data
samal <- fread("../data/samal_data.csv")
samal$time <- samal$time*24

## el hussein data
el.hussein <- fread("../data/elHussein_data.csv")
el.hussein$time <- el.hussein$time*24
el.hussein$gap <- el.hussein$gap*24

## Combine datasets
combined <- rbind(samal,el.hussein)
combined.runtimes <- combined[,.(runtime=max(time)),by=gap]

## plot elHussein data
tiff("../figures/hussein_data.tif",res=600,
     compression="lzw",height=400*9,width=400*12,pointsize=8)
par(mfrow=c(2,1),mar=c(4.1,4.1,3.1,2.1),
    oma=c(3,0,3,0))
cols=colorRampPalette(brewer.pal(n=3,name="Reds"))(10)
plot(-1,-1,xlim=c(0,16),ylim=c(0,1),xlab="Time (days)",
     ylab="Proportion reassortant",bty="n",main="Varying time of measurement")
points(samal$time/24, samal$prop_reassortant,col=cols[1],pch=1,lwd=2)
points(el.hussein[gap==0]$time/24, el.hussein[gap==0]$prop_reassortant,col=cols[1],pch=19)
points(el.hussein[gap==1*24]$time/24, el.hussein[gap==1*24]$prop_reassortant,col=cols[2],pch=19)
points(el.hussein[gap==3*24]$time/24, el.hussein[gap==3*24]$prop_reassortant,col=cols[4],pch=19)
points(el.hussein[gap==5*24]$time/24, el.hussein[gap==5*24]$prop_reassortant,col=cols[6],pch=19)
points(el.hussein[gap==7*24]$time/24, el.hussein[gap==7*24]$prop_reassortant,col=cols[8],pch=19)
points(el.hussein[gap==9*24]$time/24, el.hussein[gap==9*24]$prop_reassortant,col=cols[10],pch=19)
legend("topleft",legend=paste(c(0,3,6,9),"day gap"),bty="n",col=cols[c(0,3,6,9)+1],
       pch=19)
legend("bottomleft",legend=c("Samal et al.","El Hussein et al."),bty="n",
       col="black",
       pch=c(1,19))
cols=colorRampPalette(brewer.pal(n=3,name="Blues"))(12)
plot(-1,-1,xlim=c(0,16),ylim=c(0,1),xlab="Gap (days)",
     ylab="Proportion reassortant",bty="n",main="Varying gap between infections")
points(el.hussein[time==5*24]$gap/24, el.hussein[time==5*24]$prop_reassortant,col=cols[1],pch=19)
points(rep(0,nrow(samal[time==7*24]))/24, samal[time==7*24]$prop_reassortant,col=cols[3],pch=1,lwd=2)
points(el.hussein[time==7*24]$gap/24, el.hussein[time==7*24]$prop_reassortant,col=cols[3],pch=19)
points(el.hussein[time==9*24]$gap/24, el.hussein[time==9*24]$prop_reassortant,col=cols[5],pch=19)
points(rep(0,nrow(samal[time==10*24]))/24, samal[time==10*24]$prop_reassortant,col=cols[9],pch=1,lwd=2)
points(el.hussein[time==10*24]$gap/24, el.hussein[time==10*24]$prop_reassortant,col=cols[6],pch=19)
points(rep(0,nrow(samal[time==11*24]))/24, samal[time==11*24]$prop_reassortant,col=cols[1],pch=1,lwd=2)
points(el.hussein[time==11*24]$gap/24, el.hussein[time==11*24]$prop_reassortant,col=cols[7],pch=19)
points(rep(0,nrow(samal[time==13*24]))/24, samal[time==13*24]$prop_reassortant,col=cols[9],pch=1,lwd=2)
points(rep(0,nrow(samal[time==15*24]))/24, samal[time==15*24]$prop_reassortant,col=cols[11],pch=1,
       lwd=2)
points(rep(0,nrow(samal[time==16*24]))/24, samal[time==16*24]$prop_reassortant,col=cols[12],pch=1,
       lwd=2)
legend("topright",legend=paste(c(6,9,12,15),"days since first infection"),bty="n",
       col=cols[c(6,9,12,15)-4],
       pch=19)
legend("bottomright",legend=c("Samal et al.","El Hussein et al."),bty="n",
       col="black",
       pch=c(1,19))
##mtext(" (Days)",side=1,line=3,cex=3/4)
dev.off()

## set up parameters and state variable
state.coinf[["V.ba"]] <- floor(initial.titre*log(2) + 0.5)
state.coinf[["V.bb"]] <- 0
withlike <- parms[["withlike"]]
p.db.inc <- parms[["p.db.inc"]]
par <- c(withlike=withlike,p.db.inc=p.db.inc)

## Plot the time series
if (plot.timeseries) {
    pdf("../figures/coinfection_timeseries.pdf")
    parms.temp <- parms
    parms.temp["p.ma"] <- parms[["p.m"]]
    parms.temp["p.mb"] <- parms[["p.m"]]
    parms.temp["p.mr"] <- parms[["p.m"]]
    parms.temp["p.sa"] <- parms[["p.s"]]
    parms.temp["p.sb"] <- parms[["p.s"]]
    parms.temp["p.sr"] <- parms[["p.s"]]
    par(mfrow=c(3,2))
    for (second.intro in unique(combined.runtimes$gap)) {
        out.temp <- calc.R.variable(state.coinf,parms.temp,
                                    fitting.parms=par,runtime=max(combined.runtimes$runtime),
                                    second.intro=second.intro)
        len <- length(out.temp$V.tot)
        V.adj <- (prop.positive.vec[1:len] - prop.disseminated) * out.temp$V.tot + prop.disseminated * out.temp$V.tot.pos
        R.adj <- (prop.positive.vec[1:len] - prop.disseminated) * out.temp$R.tot + prop.disseminated * out.temp$R.tot.pos
        plot(seq(0,len-1),V.adj/log(2),type='l',log='y',col="black",lwd=2.5,
             xlab = "Time since first infection (hours)", ylab=expression(Viral~load~(TCID[50])),
             bty="n",las=1,xaxs="i",main=paste(second.intro," hour gap between infections"),
             xlim=c(0,max(combined.runtimes$runtime)))
        lines(seq(0,len-1),R.adj/log(2),col=adjustcolor("red",0.6),lwd=2.5)
        abline(v=second.intro,lty=2)
        if (second.intro == 0) {
            legend("bottomright",legend=c("Total viral load","Reassortant viral load"),
                   col=c("black",adjustcolor("red",0.6)),bty="n",lwd=2.5)
        }
    }
    dev.off()
}

## Retrieve GAM fit output
load("./gam_output_and_data.RData",verbose=TRUE)

## Get reassortment proportion with gap
if (plot.reassortment) {
    R.adjs <- list()
    V.adjs <- list()
    parms.temp <- parms
    parms.temp["p.ma"] <- parms[["p.m"]]
    parms.temp["p.mb"] <- parms[["p.m"]]
    parms.temp["p.mr"] <- parms[["p.m"]]
    parms.temp["p.sa"] <- parms[["p.s"]]
    parms.temp["p.sb"] <- parms[["p.s"]]
    parms.temp["p.sr"] <- parms[["p.s"]]
    ##out <- calc.R.full(state.coinf,parms,fitting.parms=par)
    for (second.intro in unique(combined.runtimes$gap)){ 
        out.temp <- calc.R.variable(state.coinf,parms.temp,fitting.parms=par,
                                    second.intro=second.intro,runtime=15*24)
        ##out.temp <- out[[as.character(second.intro)]]
        len <- length(out.temp$V.tot)
        V.adj <- (prop.positive.vec[1:len] - prop.disseminated) * out.temp$V.tot + prop.disseminated * out.temp$V.tot.pos
        V.adjs[[as.character(second.intro)]] <- V.adj
        R.adj <- (prop.positive.vec[1:len] - prop.disseminated) * out.temp$R.tot + prop.disseminated * out.temp$R.tot.pos
        R.adjs[[as.character(second.intro)]] <- R.adj
    }
    ## Plot it
    cols <- brewer.pal(n=6,name="Reds")
    pdf("../figures/reassortment_fitted_parms.pdf",width=7,height=3.5)
    second.intros <- unique(combined.runtimes$gap)
    par(mfrow=c(1,2),mar=c(5.1,4.1,3.1,1.1))
    ii <- 0
    for (tt in seq(5,15,2)) {
        ii <- ii + 1
        gam.pred <- predict(g.log,newdata=data.frame(time=tt,gap=seq(0,tt,1/24)),se.fit=TRUE,type='link')
        if (ii == 1) {
            plot(seq(0,tt,1/24),plogis(gam.pred$fit),type='l',
                 xlim=c(0,10),ylim=c(0,1),
                 lwd=2,xlab="Second infection day",ylab="Proportion reassortant",
                 col=cols[ii],lty=1,bty="n",las=1,yaxs="i",
                 main="Generalized additive model")
        } else {
            lines(seq(0,tt,1/24),plogis(gam.pred$fit),lwd=2,
                  col=cols[ii],lty=1)
        }
    }
    legend("topright",title="Days since 1st infection",bty="n",
           legend=seq(5,15,2),col=cols,lwd=2,lty=1)
    ii <- 0
    for (tt in seq(5,15,2)) {
        ii <- ii+1
        second.intros.temp <- second.intros[second.intros < tt *24]
        if (ii == 1) {
            plot(second.intros.temp/24,
                 sapply(R.adjs,function(x)x[[tt*24+1]])[as.character(second.intros.temp)]/
                 sapply(V.adjs,function(x)x[[tt*24+1]])[as.character(second.intros.temp)],
                 col=cols[ii],lwd=2,type='l',
                 xlim=c(0,10),ylim=c(0,1),
                 xlab="Second infection day",ylab="",
                 bty="n",las=1,yaxs="i",
                 main="Dynamical model")
        }else{
            lines(second.intros.temp/24,
                  sapply(R.adjs,function(x)x[[tt*24+1]])[as.character(second.intros.temp)]/
                  sapply(V.adjs,function(x)x[[tt*24+1]])[as.character(second.intros.temp)],
                  col=cols[ii],lwd=2)
        }
    }
    dev.off()
}
    
## Parameter sweeps
if (regenerate.parmsweep) {
    epsilon.ms <- exp(seq(log(parms["epsilon.m"]/sqrt(10)),
                          log(parms["epsilon.m"]*sqrt(10)),
                          length.out=201))
    epsilon.ss <- exp(seq(log(parms["epsilon.s"]/sqrt(10)),
                          log(parms["epsilon.s"]*sqrt(10)),
                          length.out=length(epsilon.ms)))
    p.folds <- 10^seq(log(1/10,10),log(10,10),0.05)
    runtime <- 10 * 24
    second.intros <- c(0,3 * 24)
    target.0day <- matrix(NA,nrow=length(epsilon.ms),ncol=length(p.folds))
    target.3day <- target.0day; target.ratio <- target.0day
    parms.temp <- parms
    for (ii in seq(1:length(epsilon.ms))) {
        print(ii)
        for (jj in seq(1:length(p.folds))) {
            parms.temp["epsilon.m"] <- epsilon.ms[ii]; parms.temp["epsilon.s"] <- epsilon.ss[ii]
            parms.temp["p.ma"] <- parms[["p.m"]] / sqrt(p.folds[jj])
            parms.temp["p.mb"] <- parms[["p.m"]] * sqrt(p.folds[jj])
            parms.temp["p.mr"] <- parms[["p.m"]] #max(parms.temp[["p.mb"]],parms.temp[["p.mb"]]) 
            parms.temp["p.sa"] <- parms[["p.s"]] / sqrt(p.folds[jj])
            parms.temp["p.sb"] <- parms[["p.s"]] * sqrt(p.folds[jj])
            parms.temp["p.sr"] <- parms[["p.s"]] #max(parms.temp[["p.sb"]],parms.temp[["p.sb"]])
            out.temp.1 <- calc.R.variable(state.coinf,
                                          parms.temp,fitting.parms=par,
                                          runtime=runtime,second.intro=second.intros[1])
            out.temp.2 <- calc.R.variable(state.coinf,
                                          parms.temp,fitting.parms=par,
                                          runtime=runtime,second.intro=second.intros[2])
            V.adj.1 <- (prop.positive.vec[runtime + 1] - prop.disseminated) * out.temp.1$V.tot[runtime + 1] + prop.disseminated * out.temp.1$V.tot.pos[runtime + 1]
            V.adj.2 <- (prop.positive.vec[runtime + 1] - prop.disseminated) * out.temp.2$V.tot[runtime + 1] + prop.disseminated * out.temp.2$V.tot.pos[runtime + 1]
            R.adj.1 <- (prop.positive.vec[runtime + 1] - prop.disseminated) * out.temp.1$R.tot[runtime + 1] + prop.disseminated * out.temp.1$R.tot.pos[runtime + 1]
            R.adj.2 <- (prop.positive.vec[runtime + 1] - prop.disseminated) * out.temp.2$R.tot[runtime + 1] + prop.disseminated * out.temp.2$R.tot.pos[runtime + 1]
            target.0day[ii,jj] <- R.adj.1 / V.adj.1
            target.3day[ii,jj] <- R.adj.2 / V.adj.2
            target.ratio[ii,jj] <- R.adj.1 * V.adj.2 / R.adj.2 / V.adj.1
        }
    }
    save(epsilon.ms,epsilon.ss,p.folds,
         target.0day,target.3day,target.ratio,file="sweep_epsilon_pfold_reassortment.RData")
}
 
if (plot.parmsweep) {
    load("sweep_epsilon_betafold_reassortment.RData")

    resolution <- 1200
    
    tiff("../figures/sweep_epsilon_betafold_reassortment_ratio.tif",res=resolution,compression="lzw",
         width=resolution*5.1,height=resolution*3.7)
    par(mai=c(0.22,0.62,0.42,0.42))
    nlvls <- 21
    ma <- max(abs(range(log(target.ratio,10),na.rm=TRUE)))
    lvls <- seq(-ma,ma,length.out=nlvls)
    cols.diverging <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(nlvls-1))
    filled.contour(log(epsilon.ms,10),log(beta.folds,10),
                   log(apply(target.ratio,c(1,2),function(x)max(x,0)),10),
                   ## plot.title = title(xlab="Eclipse phase of midgut cells (hours)",
                   ##                    ylab="Fold difference in viral production rates"),
                   axes=FALSE,
                   plot.axes={axis(1,at=log(c(20,50,100,200,500),10),
                                   labels=FALSE);
                                   axis(2,at=log(c(0.1,0.2,0.5,1.0,2,5,10),10),
                                        labels=FALSE);
                                   abline(h=0,lwd=2,lty=2);
                                   abline(v=log(parms["epsilon.m"],10),lwd=2,lty=3)},
                   key.axes=axis(4,at=pretty(log(apply(target.ratio,c(1,2),function(x)max(x,0)),10)),
                                 labels=round(10^pretty(log(apply(target.ratio,c(1,2),
                                                                  function(x)max(x,0)),10)),5)),
                   col=cols.diverging,
                   key.title={par(cex.main=0.8);title("Co-infection\n/ sequential")},
                   levels=lvls)
    ## mtext(side = 3, line = 1, adj = -0.05, 'A', font = 2)
    dev.off()

    tiff("../figures/sweep_slices_reassortment_eclipse.tif",res=resolution,compression="lzw",
         width=resolution*3.9,height=resolution*1.9)
    par(mar=c(4.1,3.1,2.1,0.4))
    ## par(mfrow=c(2,1))
    plot(epsilon.ms,target.ratio[,which(beta.folds==1)],bty="n",
         xlab="Eclipse phase of midgut cells (hours)",
         type='l',lwd=3,
         ylab="",las=1,log="xy",xaxs="i",xaxt="n")
    axis(1,at=c(20,50,100,200,500))
    abline(h=1,lty=4,lwd=2)
    abline(v=parms["epsilon.m"],lty=3,lwd=2)
    ## mtext(side = 3, line = 1, adj = -0.14, 'B', font = 2)
    dev.off()

    tiff("../figures/sweep_slices_reassortment_beta.tif",res=resolution,compression="lzw",
         width=resolution*1.9,height=resolution*4)
    par(mar=c(2.2,4.1,1.75,1.1))
    plot(target.ratio[which.min(abs(epsilon.ms - parms["epsilon.m"])),],
         beta.folds,
         bty="n",ylab=expression(paste("Fold difference in ",beta)),
         type='l',lwd=3,
         xlab="",las=1,log="xy",yaxs="i")
    abline(v=1,lty=4,lwd=2)
    abline(h=1,lty=2,lwd=2)
    mtext(side = 3, line = 0.8, adj = -0.6, 'C', font = 2)
    dev.off()

    pdf("../figures/combined_sweep_ratio_beta.pdf",width=7,height=5.5)
    layout(matrix(c(rep(rev(c(rep(1,52),rep(3,18))),37),
                    rep(c(rep(3,18),rep(2,39),rep(5,13)),3),
                    rep(c(rep(4,18),rep(2,39),rep(5,13)),15)),
                  byrow=TRUE,ncol=70))
    par(mar = c(0,0,0,0))
    im <- load.image("../figures/sweep_epsilon_betafold_reassortment_ratio.tif")
    plot(im,axes=FALSE)
    ##par(mar = c(0,0.1,1.3,0.1))
    im <- load.image("../figures/sweep_slices_reassortment_eclipse.tif")
    plot(im,axes=FALSE)
    im <- load.image("../figures/sweep_slices_reassortment_beta.tif")
    plot(im,axes=FALSE)
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
    text(0.2,0.6,
         "Ratio of the\nprop. reassortant\nfollowing\nco-infection\nto sequential\ninfections\nseparated by 3 days",
         cex=4/3,adj=0)
    dev.off()

    tiff("../figures/sweep_epsilon_betafold_reassortment_0day.tif",res=resolution,compression="lzw",
         width=resolution*5.1,height=resolution*3.7)
    par(mai=c(0.22,0.62,0.42,0.42))
    filled.contour(log(epsilon.ms,10),log(beta.folds,10),target.0day,
                   ## plot.title = title(xlab="Eclipse phase of midgut cells (hours)",
                   ##                    ylab="Fold difference in viral production rates"),
                   axes=FALSE,
                   plot.axes={axis(1,at=log(c(20,50,100,200,500),10),
                                   labels=FALSE);
                                   axis(2,at=log(c(0.1,0.2,0.5,1.0,2,5,10),10),
                                        labels=FALSE);
                                   abline(h=0,lwd=2,lty=2);
                                   abline(v=log(parms["epsilon.m"],10),lwd=2,lty=3)},
                   key.axes=axis(4,at=pretty(target.0day)),
                   key.title={par(cex.main=0.8);title("Proportion\nreassortant")},
                   color=colorRampPalette(brewer.pal(9,"Reds")))
    ## mtext(side = 3, line = 1, adj = -0.05, 'A', font = 2)
    dev.off()

    tiff("../figures/sweep_slices_0day_eclipse.tif",res=resolution,compression="lzw",
         width=resolution*3.9,height=resolution*1.9)
    par(mar=c(4.1,3.1,2.1,0.4))
    ## par(mfrow=c(2,1))
    plot(epsilon.ms,target.0day[,which(beta.folds==1)],bty="n",
         xlab="Eclipse phase of midgut cells (hours)",
         type='l',lwd=3,
         ylab="",las=1,log="x",xaxs="i",xaxt="n")
    axis(1,at=c(20,50,100,200,500))
    abline(v=parms["epsilon.m"],lty=3,lwd=2)
    ## mtext(side = 3, line = 1, adj = -0.14, 'B', font = 2)
    dev.off()

    tiff("../figures/sweep_slices_0day_beta.tif",res=resolution,compression="lzw",
         width=resolution*1.9,height=resolution*4)
    par(mar=c(2.2,4.1,1.75,1.1))
    plot(target.0day[which.min(abs(epsilon.ms - parms["epsilon.m"])),],
         beta.folds,
         bty="n",ylab=expression(paste("Fold difference in ",beta)),
         type='l',lwd=3,
         xlab="",las=1,log="xy",yaxs="i",xaxt="n")
    abline(h=1,lty=2,lwd=2)
    axis(1,at=10^(-10:0))
    ## mtext(side = 3, line = 0.8, adj = -0.6, 'C', font = 2)
    dev.off()
    
    pdf("../figures/combined_sweep_0day_beta.pdf",width=7,height=5.5)
    layout(matrix(c(rep(rev(c(rep(1,52),rep(3,18))),37),
                    rep(c(rep(3,18),rep(2,39),rep(5,13)),3),
                    rep(c(rep(4,18),rep(2,39),rep(5,13)),15)),
                  byrow=TRUE,ncol=70))
    par(mar = c(0,0,0,0))
    im <- load.image("../figures/sweep_epsilon_betafold_reassortment_0day.tif")
    plot(im,axes=FALSE)
    ## mtext(side = 3, line = 0, adj = 0.04, 'A', font = 2)
    ##par(mar = c(0,0.1,1.3,0.1))
    im <- load.image("../figures/sweep_slices_0day_eclipse.tif")
    plot(im,axes=FALSE)
    im <- load.image("../figures/sweep_slices_0day_beta.tif")
    plot(im,axes=FALSE)
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
    text(0.2,0.6,
         "Proportion\nreassortant\non day 10\nfollowing\nco-infection",
         cex=4/3,adj=0)
    dev.off()

    tiff("../figures/sweep_epsilon_betafold_reassortment_3day.tif",res=resolution,compression="lzw",
         width=resolution*5.1,height=resolution*3.7)
    par(mai=c(0.22,0.62,0.42,0.42))
    filled.contour(log(epsilon.ms,10),log(beta.folds,10),target.3day,
                   ## plot.title = title(xlab="Eclipse phase of midgut cells (hours)",
                   ##                    ylab="Fold difference in viral production rates"),
                   axes=FALSE,
                   plot.axes={axis(1,at=log(c(10,20,50,100,200,500),10),
                                   labels=FALSE);
                                   axis(2,at=log(c(0.1,0.2,0.5,1.0,2,5,10),10),
                                        labels=FALSE);
                                   abline(h=0,lwd=2,lty=2);
                                   abline(v=log(parms["epsilon.m"],10),lwd=2,lty=3)},
                   key.axes=axis(4,at=pretty(target.3day)),
                   key.title={par(cex.main=0.8);title("Proportion\nreassortant")},
                   color=colorRampPalette(brewer.pal(9,"Reds")))
    ## mtext(side = 3, line = 1, adj = -0.05, 'A', font = 2)
    dev.off()

    tiff("../figures/sweep_slices_3day_eclipse.tif",res=resolution,compression="lzw",
         width=resolution*3.9,height=resolution*1.9)
    par(mar=c(4.1,3.1,2.1,0.4))
    ## par(mfrow=c(2,1))
    plot(epsilon.ms,target.3day[,which(beta.folds==1)],bty="n",
         xlab="Eclipse phase of midgut cells (hours)",
         type='l',lwd=3,
         ylab="",las=1,log="x",xaxs="i",xaxt="n")
    abline(v=parms["epsilon.m"],lty=3,lwd=2)
    axis(1,at=c(10,20,50,100,200,500))
    ## mtext(side = 3, line = 1, adj = -0.14, 'B', font = 2)
    dev.off()

    tiff("../figures/sweep_slices_3day_beta.tif",res=resolution,compression="lzw",
         width=resolution*1.9,height=resolution*4)
    par(mar=c(2.2,4.1,1.75,1.1))
    plot(target.3day[which.min(abs(epsilon.ms - parms["epsilon.m"])),],
         beta.folds,
         bty="n",ylab=expression(paste("Fold difference in ",beta)),
         type='l',lwd=3,
         xlab="",las=1,log="xy",yaxs="i")
    abline(h=1,lty=2,lwd=2)
    ##mtext(side = 3, line = 0.8, adj = -0.6, 'C', font = 2)
    dev.off()
    
    pdf("../figures/combined_sweep_3day_beta.pdf",width=7,height=5.5)
    layout(matrix(c(rep(rev(c(rep(1,52),rep(3,18))),37),
                    rep(c(rep(3,18),rep(2,39),rep(5,13)),3),
                    rep(c(rep(4,18),rep(2,39),rep(5,13)),15)),
                  byrow=TRUE,ncol=70))
    par(mar = c(0,0,0,0))
    im <- load.image("../figures/sweep_epsilon_betafold_reassortment_3day.tif")
    plot(im,axes=FALSE)
    ## mtext(side = 3, line = 0, adj = 0.04, 'A', font = 2)
    ##par(mar = c(0,0.1,1.3,0.1))
    im <- load.image("../figures/sweep_slices_3day_eclipse.tif")
    plot(im,axes=FALSE)
    im <- load.image("../figures/sweep_slices_3day_beta.tif")
    plot(im,axes=FALSE)
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
    text(0.2,0.6,
         "Proportion\nreassortant\non day 10\nfollowing\nsequential infections\nseparated\nby 3 days",
         cex=4/3,adj=0)
    dev.off()

}

## Sobol analysis
if (plot.sobol){
    load("./sobol_sweep_combined_output.RData",verbose=TRUE)
    indices.0day <- sobol_indices(Y=output.combined[,16],N=1e5,params=colnames(output.combined)[1:15])
    indices.3day <- sobol_indices(Y=output.combined[,17],N=1e5,params=colnames(output.combined)[1:15])
    indices.ratio <- sobol_indices(Y=output.combined[,18],N=1e5,params=colnames(output.combined)[1:15])

    pdf("../figures/sobol_pies.pdf",width=10,height=4,pointsize=14)
    layout(t(c(1,1,2,2,3,3)))
    par(mar=c(0,0,3,0))
    pie(c(pmax(0,unlist(indices.0day$results[1:15,1])),
          1-pmax(0,sum(unlist(indices.0day$results[1:15,1])))),
        labels=c(rep("",4),expression(mu[m]),rep("",6),expression(mu[s]),
                 "",expression(omega),"","Interactions"),
        col=viridis(16),
        main="Simultaneous co-infection")
    plot(-1,-1,bty="n",ylim=c(0,1),xlim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="")
    legend("center",legend=c(expression(beta[m] * T[m]),expression(d[m]),
                             expression(p[m]),expression(T["0,m"]),expression(mu[m]),
                             expression(epsilon[m]),expression(c[s]),
                             expression(beta[s]),expression(d[s]),expression(p[s]),
                             expression(T["0,s"]),
                             expression(mu[s]),expression(epsilon[s]),
                             expression(omega),expression(P[db]^(inc)),
                             "Interactions"),
           fill=viridis(16),bty="n",ncol=2)
    pie(c(pmax(0,unlist(indices.3day$results[1:15,1])),
              1-pmax(0,sum(unlist(indices.3day$results[1:15,1])))),
        labels=c(rep("",4),expression(mu[m]),"",expression(c[s]),
                 rep("",4),expression(mu[s]),"",expression(omega),"","Interactions"),
        col=viridis(16),main="Second infection on day 3")
    dev.off()
    pdf("../figures/sobol_pies_higher.pdf",width=10,height=4,pointsize=14)
    layout(t(c(1,1,2,2,3,3)))
    par(mar=c(0,0,3,0))
    pie(pmax(0,unlist(indices.0day$results[16:30,1])),
        labels=c(rep("",4),expression(mu[m]),rep("",6),expression(mu[s]),"",expression(omega),""),
        col=viridis(16)[1:15],
        main="Simultaneous co-infection")
    plot(-1,-1,bty="n",ylim=c(0,1),xlim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="")
    legend("center",legend=c(expression(beta[m] * T[m]),expression(d[m]),
                             expression(p[m]),expression(T["0,m"]),expression(mu[m]),
                             expression(epsilon[m]),expression(c[s]),
                             expression(beta[s]),expression(d[s]),expression(p[s]),
                             expression(T["0,s"]),
                             expression(mu[s]),expression(epsilon[s]),
                             expression(omega),expression(P[db]^(inc))),
           fill=viridis(16)[1:15],bty="n",ncol=2)
    pie(pmax(0,unlist(indices.3day$results[16:30,1])),
        labels=c(rep("",4),expression(mu[m]),"",expression(c[s]),
                 rep("",4),expression(mu[s]),"",expression(omega),""),
        col=viridis(16)[1:15],main="Second infection on day 3")
    dev.off()

    second.indices.0day <- sobol_indices(Y=output.combined[,16],N=1e5,
                                         params=colnames(output.combined)[1:15],order="second")
    second.indices.3day <- sobol_indices(Y=output.combined[,17],N=1e5,
                                         params=colnames(output.combined)[1:15],order="second")
}
