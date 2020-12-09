### General functions
temperature.curve <- function(time,mean,amp,period=1) mean + amp*sin(time*2*pi/period)
integrate <- function(integrand, Dt) {
    Dt*(sum(integrand) - (integrand[1] + integrand[length(integrand)])/2)
}
rate.prob <- function(rate.fn, Tvec, Dt) {
    rates <- rate.fn(Tvec)
    integral <- integrate(rates, Dt)
    return(1-exp(-integral))
}
period.prob <- function(period.fn, Tvec, Dt) {
    rates <- 1/period.fn(Tvec)
    integral <- integrate(rates,Dt)
    return(1-exp(-integral))
}
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    ## dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(4, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}
growth.rate <- function(period.function, threshold, T) {
    log(threshold)/period.function(T)
}
