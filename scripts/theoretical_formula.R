## Mortality:
##    https://academic.oup.com/jme/article-abstract/37/5/675/954678?redirectedFrom=fulltext (Fig. 9)
## Incubation: ditto, figure 2
## Alternative mortality:
##    https://onlinelibrary.wiley.com/doi/epdf/10.1046/j.1365-2915.2002.00357.x (Table 3)
## Alternative incubation period estimates:
##    https://onlinelibrary.wiley.com/doi/epdf/10.1046/j.1365-2915.2002.00357.x (Fig. 2)
##    https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2915.1995.tb00119.x 
## Biting rates: https://pubmed.ncbi.nlm.nih.gov/19496435/
##    Or try formula in https://onlinelibrary.wiley.com/doi/full/10.1046/j.1365-2915.2002.00357.x
##    which comes from Mullens & Holbrook
## Vector competence:
##    https://onlinelibrary.wiley.com/doi/epdf/10.1046/j.1365-2915.2002.00357.x (Fig. 3)
##    Not T-dependent, around 9.76% for BTV10 and 12.4% for BTV16
## all time units are days
## Use the Wittmann paper which gathers estimates of all quantities in one place.
library(viridis)

recalculate <- FALSE

### Parameters - assumed
Dt <- 1/24
tmax <- 61
tvec <- seq(0,tmax,Dt)
Tmean <- 21
Tamp <- 11

###Parameters - literature derived
vector.competence.10 <- sin(18.2343*pi/180)^2
vector.competence.16 <- sin(20.6322*pi/180)^2

###Load functions
source("./general_functions.R")
source("./btv_culicoides_functions.R")

## Plot different incubation periods
Tvec.inc <- seq(20,30,0.1)
tiff("../figures/incubation_comparison.tif",res=600,compression="lzw",height=600*6,width=600*6)
plot(Tvec.inc,incubation.period(Tvec.inc),type='l',
     ylim=range(c(incubation.period(Tvec.inc),
                  incubation.period.10(Tvec.inc),
                  incubation.period.16(Tvec.inc))),
     ylab="Incubation period (days)",
     xlab="Temperature (C)")
lines(Tvec.inc,incubation.period.10(Tvec.inc),col="red")
lines(Tvec.inc,incubation.period.16(Tvec.inc),col="gray")
abline(h=c(11,14),col="light gray",lty="dotted")
abline(v=c(20,25,30),col="light gray",lty="dashed")
legend(24.8,16.4,legend=c("Gerry & Mullens","Wittmann et al. BTV10","Wittmann et al. BTV16"),
       col=c("black","red","gray"),lty="solid",bty="n")
dev.off()

## 1. and 2. Plot probabilities/rates
Tvec <- temperature.curve(tvec,Tmean,Tamp)
survival.probs <- vector(mode="numeric",length=length(Tvec))
incubation.probs <- vector(mode="numeric",length=length(Tvec))
reassortment.probs <- incubation.probs

for (i in 1:length(Tvec)) {
    print(i)
    survival.probs[i] <- 1-rate.prob(hazard.lowRH,Tvec[1:i],Dt)
    incubation.probs[i] <- period.prob(incubation.period.10,Tvec[1:i],Dt)
    ## placeholder for reassortment - assume rate 1/10th of incubation rate
    reassortment.probs[i] <- rate.prob(function(x){0.1/incubation.period(x)},Tvec[1:i],Dt)
}

tiff("../figures/rates_and_probabilities.tif",res=600,compression="lzw",height=600*6,width=600*6)
par(mfrow=c(5,1),mar=c(2.1,4.1,1.1,4.1),oma=c(3,0,3,0))
plot(tvec, survival.probs, type='l', ylim=c(0,1),xlab="",ylab="Survival prob.")
par(new=TRUE)
plot(tvec,hazard.lowRH(Tvec),type='l',col="gray",yaxt="n",ylab="",xaxt="n",xlab="")
axis(4,labels=T,col="gray",col.ticks="gray",col.axis="gray")
mtext("Mortality rate",4,3,col="gray",cex=2/3)

plot(tvec, incubation.probs, type='l', ylim=c(0,1),xlab="",ylab="Incubation prob.")
par(new=TRUE)
plot(tvec,incubation.period.10(Tvec),type='l',col="gray",yaxt="n",ylab="",xaxt="n",xlab="")
axis(4,labels=T,col="gray",col.ticks="gray",col.axis="gray")
mtext("Incubation period",4,3,col="gray",cex=2/3)

plot(tvec, reassortment.probs, type='l', ylim=c(0,1),xlab="",ylab="Reassortment prob.")
par(new=TRUE)
plot(tvec,0.1/incubation.period(Tvec),type='l',col="gray",yaxt="n",ylab="",xaxt="n",xlab="")
axis(4,labels=T,col="gray",col.ticks="gray",col.axis="gray")
mtext("Reassortment rate",4,3,col="gray",cex=2/3)

plot(tvec, biting.rate(Tvec), type='l',xlab="",ylab="Expected bites / day")

reassortant.bites <- biting.rate(Tvec)*survival.probs*incubation.probs*reassortment.probs*vector.competence.10
plot(tvec, reassortant.bites,
     type='l',xlab="",ylab="Infectious bites")
mtext("with reassortant / day",2,2,cex=2/3)
polygon(c(tvec,rev(tvec)),c(reassortant.bites,rep(0,length(reassortant.bites))),col="red",border=F)
dev.off()

## 3. Sum reassortant.bites
print(integrate(reassortant.bites,Dt))

total.reassortant.bites <- function(tmax,Dt,Tmean,Tamp,ratio) {
    tvec <- seq(0,tmax,Dt)
    Tvec <- temperature.curve(tvec,Tmean,Tamp)
    survival.probs <- vector(mode="numeric",length=length(Tvec))
    incubation.probs <- vector(mode="numeric",length=length(Tvec))
    reassortment.probs <- incubation.probs
    for (i in 1:length(Tvec)) {
        survival.probs[i] <- 1-rate.prob(hazard.lowRH,Tvec[1:i],Dt)
        incubation.probs[i] <- period.prob(incubation.period.10,Tvec[1:i],Dt)
        ## placeholder for reassortment - assume rate 1/10th of incubation rate
        reassortment.probs[i] <- rate.prob(function(x){ratio/incubation.period.10(x)},Tvec[1:i],Dt)
    }
    reassortant.bites <- biting.rate(Tvec)*survival.probs*incubation.probs*reassortment.probs*vector.competence.10
    return(integrate(reassortant.bites,Dt))
}

###contour plot
Tmeans <- 10:30
Tamps <- 0:15
Ebites <- array(NA,dim=c(length(Tmeans),length(Tamps)))
tmax <- 100
Dt <- 1/24
ratio <- 0.1

if (recalculate) {
    for (mm in 1:length(Tmeans)) {
        print(Tmeans[mm])
        for (aa in 1:length(Tamps)) {
            Ebites[mm,aa] <- total.reassortant.bites(100,Dt,Tmeans[mm],Tamps[aa],ratio)
        }
    }
} else {
    load("../output/Ebites.RData")
}

## image(Tmeans,Tamps,Ebites,col=magma(100),xlab="Mean (C)",ylab="Amplitude (C)")
tiff("../figures/temperature_contour.tif",res=600,compression="lzw",height=600*6,width=600*6)
par(oma=c(0,0,0,1))
filled.contour(Tmeans,Tamps,Ebites,color.palette=magma,
               xlab=expression(paste("Mean temperature("^degree,"C)")),
               ylab=expression(paste("Diurnal temperture range ("^degree,"C)")))
mtext("Expected reassortant bites",side=4,line=2)
dev.off()

## different ratios of reassortment
ratios <- 10^seq(-3,0,0.2)
Tmeans <- 10:30
Tamp <- 10
if (recalculate) {
    Ebites.rat <- array(NA,c(length(ratios),length(Tmeans)))
    for (rr in 1:length(ratios)) {
        print(ratios[rr])
        for (mm in 1:length(Tmeans)) {
            Ebites.rat[rr,mm] <- total.reassortant.bites(100,Dt,Tmeans[mm],Tamp,ratios[rr])
        }
    }
}

tiff("../figures/reassortment_line.tif",res=600,compression="lzw",height=600*6,width=600*6)
num.rows=6
par(oma=c(0,0,0,2))
layout(matrix(c(rep(1,num.rows*(num.rows-1)), rep(2,,num.rows)), nrow = num.rows))
cols <- plasma(ncol(Ebites.rat))
plot(ratios,Ebites.rat[,1],type='l',lwd=2,col=cols[1],las=1,xaxs="i",yaxs="i",
     ylim=c(0,max(Ebites.rat,na.rm=T)),xlab="Ratio reassortment:incubation rate",
     ylab="Expected reassortant bites",bty="n",las=1,xaxs="i",yaxs="i")
for (ii in 2:ncol(Ebites.rat)) {
    lines(ratios,Ebites.rat[,ii],col=cols[ii],lwd=2)
}
color.bar(cols,min=min(Tmeans),max=max(Tmeans))
mtext(expression(paste("Mean temperature("^degree,"C)")),side=4,line=3,cex=2/3)
dev.off()

tiff(file="../figures/reassortment_contour.tif",res=600,compression="lzw",height=600*6,width=600*6)
par(oma=c(0,0,0,1))
filled.contour(Tmeans,ratios,t(Ebites.rat),color.palette=magma,
               xlab=expression(paste("Mean temperature("^degree,"C)")),
               ylab="Ratio reassortment:incubation rate")
mtext("Expected reassortant bites",side=4,line=2)
dev.off()

if (recalculate) {
    save(Ebites,Ebites.rat,file="../output/Ebites.RData")
}


## reassortment based on relative growth rate
virus.threshold <- 100 #virus abundance is all relative so numbers don't really matter
initial.virus.10 <- 1
initial.virus.16 <- 1
viral.abundance.10 <- vector(mode="numeric",length=length(Tvec))
viral.abundance.16 <- viral.abundance.10
for (i in 1:length(Tvec)) {
    print(i)
    incubation.rates.10 <- 1/incubation.period.10(Tvec[1:i])
    incubation.rates.16 <- 1/incubation.period.16(Tvec[1:i])
    viral.abundance.10[i] <- initial.virus.10*virus.threshold^integrate(incubation.rates.10,Dt)
    viral.abundance.16[i] <- initial.virus.16*virus.threshold^integrate(incubation.rates.16,Dt)
}

plot(tvec,viral.abundance.16,type='l',log='y',xlab="Time (days)",
     ylab="Relative virus population",lwd=2)
lines(tvec,viral.abundance.10,lwd=2)
abline(h=virus.threshold,lty="dashed")
baseline.reassortment.rate <- 0.1
reassortants <- baseline.reassortment.rate*viral.abundance.10*viral.abundance.16
par(new=TRUE)
plot(tvec,reassortants, type='l',yaxt='n',ylab='',log='y',col='red',lty='dotted',lwd=2,
     xlab="")
