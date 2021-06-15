## set working directory
setwd('~/Documents/bluetongue_project/aim2/scripts/')

## clear existing workspacer
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(BayesianTools,psych,bbmle,deSolve,grDevices,data.table,plotrix)


## Fu data
fu <- read.csv("../data/Fu_data.csv",header=T)
fu$log10.mean.titre <- log(log(2) * 10^fu$log10.mean.titre,10) # convert to pfu
fu.hours <- fu$Time.pi..hours.

## Hemocoel only
fu.intrathoracic <- read.csv("../data/Fu_data_intrathoracic.csv",header=T)
fu.intrathoracic$titre <- fu.intrathoracic$titre * log(2) #convert to pfu
fu.intrathoracic.hours <- pmax(floor(fu.intrathoracic$day*24 + 0.5),0)

## Samal data
samal <- fread("../data/samal_data.csv")
samal$time <- samal$time*24

## Samal isolation data
samal.isolation <- fread("../data/samal_data_isolation.csv")
samal.isolation$time <- samal.isolation$time*24

## el hussein data
el.hussein <- fread("../data/elHussein_data.csv")
el.hussein$time <- el.hussein$time*24
el.hussein$gap <- el.hussein$gap*24

plot(samal$titre,samal$prop_reassortant,log="x")
cor(samal$titre,samal$prop_reassortant )
plot(samal$time,samal$titre,log="y")

plot(fu$Time.pi..hours., fu$Detection.rate,type='l',ylim=c(0,1),xlim=c(0,16*24))
plotCI(samal.isolation$time,samal.isolation$pos/samal.isolation$n,
       ui=samal.isolation$pos/samal.isolation$n + 1.96*sqrt(samal.isolation$pos*(samal.isolation$n-samal.isolation$pos)/samal.isolation$n^3),
       li=samal.isolation$pos/samal.isolation$n - 1.96*sqrt(samal.isolation$pos*(samal.isolation$n-samal.isolation$pos)/samal.isolation$n^3),add=TRUE)
