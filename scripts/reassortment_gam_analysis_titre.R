if(!require(pacman)) {install.packages("pacman");library(pacman)}
p_load(data.table,mgcv,scam)

combined <- fread("../data/elHussein_samal_data_combined.csv")
combined=combined[titre!=1]

## Fit GAm
g <- gam(log(titre,10)~s(time,k=4)+s(gap,k=4)+experiment,data=combined,method="REML")
coef(g)
summary(g)
pdf("../figures/gam_output_titre.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(g,shift=coef(g)[1],seWithMean=TRUE,
     yaxt="n",select=1,xlab="Time of measurement (days)",ylab="Plaque-forming units / fly",
     shade=TRUE,rug=FALSE,ylim=log(range(combined$titre),10))
points(log(titre,10)~time,data=combined)
par(new=TRUE)
axis(2,at=1:5,labels=10^(1:5),las=1)
plot(g,shift=coef(g)[1],seWithMean=TRUE,
     yaxt="n",select=2,xlab="Gap between infections (days)",ylab="",
     shade=TRUE,rug=FALSE,ylim=log(range(combined$titre),10))
points(log(titre,10)~gap,data=combined)
axis(2,at=1:5,labels=10^(1:5),las=1)
dev.off()

gam.check(g)

## Fit GAM without anomalous points
g.anom <- gam(titre~s(time,k=5)+s(gap,k=5)+experiment,data=combined[titre<2e5 & titre>2],method="REML")
coef(g.anom)
summary(g.anom)
plot(g.anom,pages=1,pch=1,residuals=TRUE,shift=coef(g.anom)[1])

gam.check(g.anom)

## Fit GAM on log-scale
g.log <- gam(log(titre)~s(time,k=5)+s(gap,k=5)+experiment,data=combined,method="REML")
coef(g.log)
summary(g.log)
plot(g.log,pages=1,pch=1,residuals=TRUE,shift=coef(g.log)[1])

gam.check(g.log)
## Fit GAM without anomalous points on log-scale
g.anomlog <- gam(log(titre)~s(time,k=5)+s(gap,k=5)+experiment,data=combined[titre<2e5 & titre>2],method="REML")
coef(g.anomlog)
summary(g.anomlog)
plot(g.anomlog,pages=1,pch=1,residuals=TRUE,shift=coef(g.anomlog)[1])

gam.check(g.anomlog)

## Try a log-link
## Fit GAm
g <- gam(titre~s(time,k=5)+s(gap,k=5)+experiment,data=combined,method="REML")
coef(g)
summary(g)
plot(g,pages=1,pch=1,residuals=TRUE,shift=coef(g)[1],
     seWithMean=TRUE,shade=TRUE,rug=FALSE)

gam.check(g)
