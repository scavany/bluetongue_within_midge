library(data.table)

fu.data <- fread("../data/Fu_data.csv")
names(fu.data) <- c("time","midges","log10titre","lower","upper","detection")
plot(fu.data$time,fu.data$detection,type="l")
t0 <- 30
k <- -1
l <- 0.35
lines(fu.data$time,(1-l)/(1+exp(-k*(fu.data$time - t0))) + l,col="red")
