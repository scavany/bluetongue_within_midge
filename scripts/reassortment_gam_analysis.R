if(!require(pacman)) {install.packages("pacman");library(pacman)}
p_load(data.table,mgcv,scam)

## Samal data
samal <- fread("../data/samal_data.csv")
samal$time <- samal$time*24
## samal$gap <- 0
## samal[,total:=n_plaques]

## el hussein data
el.hussein <- fread("../data/elHussein_data.csv")
el.hussein$time <- el.hussein$time*24
el.hussein$gap <- el.hussein$gap*24

## Combine data sets
combined <- rbind(samal[,.(time,gap,prop_reassortant,titre,total)],
                  el.hussein[,.(time,gap,prop_reassortant,titre,total)])
combined[,logit.reassortant := -log(1/prop_reassortant - 1)]

## Run GAM
g = gam(prop_reassortant~s(time,k=3)+s(gap,k=5),data=combined)
coef(g)
summary(g)
plot(g,residuals=TRUE,pch=1,pages=1,shift = coef(g)[1])
gam.check(g)

g.parmgap = gam(prop_reassortant~s(time,k=3)+gap,data=combined)
summary(g.parmgap)
plot(g.parmgap,residuals=TRUE,pch=1,pages=1)

##Try titre and log.titre
g.titre <- gam(titre~s(time,k=3)+s(gap,k=3),data=combined)
plot(g.titre,residuals=TRUE,pch=1,pages=1)
gam.check(g.titre)

## Try SCAM
s = scam(prop_reassortant~s(time,bs="mpi")+s(gap,bs="mpd"),data=combined)
plot(s,residuals=TRUE,pch=1,pages=1,shift = coef(s)[1])
summary(s)

## Try logistic regression
## Create data
time.temp <- rep(combined$time,combined$total)
gap.temp <- rep(combined$gap,combined$total)
n.reassort.temp <- floor(combined$prop_reassortant*combined$total+0.5)
reassort.temp <- c()
for (ii  in 1:length(n.reassort.temp)){
     reassort.temp <- c(reassort.temp,
                        rep(c(1,0),
                            c(n.reassort.temp[ii],combined$total[ii]-n.reassort.temp[ii])))
}
combined.binary <- data.frame(list(time=time.temp,gap=gap.temp,
                                   reassortant=as.logical(reassort.temp)))
setDT(combined.binary)
rm(time.temp,gap.temp,n.reassort.temp,reassort.temp)

## Fit model
model <- glm(reassortant~time+gap,family=binomial(link=logit),data=combined.binary)
summary(model)

## Do a gam on this data
combined.binary[,':='(time=time/24,gap=gap/24)]
g.log <- gam(reassortant~s(time,k=5)+s(gap,k=4),data=combined.binary,
             family=binomial,
             method = "REML")
summary(g.log)

pdf("../figures/gam_output.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(g.log,trans=plogis,shift=coef(g.log)[1],
     seWithMean=TRUE,shade=TRUE,rug=FALSE,select=1,
     xlab="Time of measurement (days)",
     ylab="Proportion reassortant")
for (t in unique(combined.binary$time)) {
    points(t,sum(combined.binary[time==t]$reassortant) /
                length(combined.binary[time==t]$reassortant))
}
plot(g.log,trans=plogis,shift=coef(g.log)[1],
     seWithMean=TRUE,shade=TRUE,rug=FALSE,select=2,
     xlab="Gap between infections (days)",
     ylab="")
for (t in unique(combined.binary$gap)) {
    points(t,sum(combined.binary[gap==t]$reassortant) /
                length(combined.binary[gap==t]$reassortant))
}
dev.off()

gam.check(g.log)

## Do a SCAM
s.log <- scam(reassortant~s(time,bs="mpi")+s(gap,bs="mpd"),
              data=combined.binary,family=binomial(link=logit))
summary(s.log)
plot(s.log,pages=1,residuals=TRUE,pch=1,trans=plogis,shift=coef(g.log)[1],
     seWithMean=TRUE,rug=FALSE,shade=TRUE)


## Boneyard
## Now fit on el Hussein and test on Samal data
time.temp <- rep(el.hussein$time,el.hussein$total)
gap.temp <- rep(el.hussein$gap,el.hussein$total)
n.reassort.temp <- floor(el.hussein$prop_reassortant*el.hussein$total+0.5)
reassort.temp <- c()
for (ii  in 1:length(n.reassort.temp)){
     reassort.temp <- c(reassort.temp,
                        rep(c(1,0),
                            c(n.reassort.temp[ii],el.hussein$total[ii]-n.reassort.temp[ii])))
}
el.hussein.binary <- data.frame(list(time=time.temp,gap=gap.temp,
                                   reassortant=as.logical(reassort.temp)))
setDT(el.hussein.binary)
rm(time.temp,gap.temp,n.reassort.temp,reassort.temp)
time.temp <- rep(samal$time,samal$total)
gap.temp <- rep(samal$gap,samal$total)
n.reassort.temp <- floor(samal$prop_reassortant*samal$total+0.5)
reassort.temp <- c()
for (ii  in 1:length(n.reassort.temp)){
     reassort.temp <- c(reassort.temp,
                        rep(c(1,0),
                            c(n.reassort.temp[ii],samal$total[ii]-n.reassort.temp[ii])))
}
samal.binary <- data.frame(list(time=time.temp,gap=gap.temp,
                                   reassortant=as.logical(reassort.temp)))
setDT(samal.binary)
rm(time.temp,gap.temp,n.reassort.temp,reassort.temp)

model <- glm(reassortant~time+gap,family=binomial(link=logit),data=el.hussein.binary)
summary(model)
predict(model,samal.binary,type="response")
sum(samal.binary$reassortant)/length(samal.binary$reassortant)
mean(predict(model,samal.binary,type="response"))
