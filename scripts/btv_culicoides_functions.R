### BTV/Culicoides specific functions
## life.expectancy <- function(T) 111.84*exp(-0.1547*T)
hazard.lowRH <- function(T) pmax(0,0.0036*T - 0.0165)
hazard.highRH <- function(T) pmax(0,0.0064*T - 0.0681)
incubation.period <- function(T) pmax(0,-1.03*T+36.79)
incubation.period.10 <- function(T) 1/(pmax(0,0.0069*T-0.0636))
incubation.period.16 <- function(T) 1/(pmax(0,0.0113*T-0.1419))
biting.rate <- function(T) 1/pmax(0,-1.98 + 0.07217*T + 2516.65/T^2)
