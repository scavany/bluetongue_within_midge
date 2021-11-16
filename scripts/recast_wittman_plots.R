## Fig. 2
pdf("../figures/wittmann_fig2_recast.pdf",width=2.8,height=4.4)
par(mfrow=c(2,1),mar=c(2.1,2.1,1.1,1.1))
Tmin <- 9.2
x <- seq(Tmin+0.1,30,0.1)
m <- 0.0069
c <- -0.0636
y <- m * x + c
plot(x,1/y,type='l',lwd=2.5,bty="n",las=1,ylim=c(0,61),xlim=c(5,30),xaxs="i",yaxs="i",
     xlab="",ylab="")
abline(v=Tmin,lty=2)
Tmin <- 12.6
x <- seq(Tmin+0.1,30,0.1)
m <- 0.0113
c <- -0.1419
y <- m * x + c
plot(x,1/y,type='l',lwd=2.5,bty="n",las=1,ylim=c(0,61),xlim=c(5,30),xaxs="i",yaxs="i",
     xlab="",ylab="")
abline(v=Tmin,lty=2)
dev.off()

