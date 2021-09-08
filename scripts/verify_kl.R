# TO DO: Explain here what this script does, and how to use it.

# Script settings.
n  <- 75
s  <- 3
b  <- -1
s0 <- c(0.01,0.1,1)
w0 <- c(0.9,0.05,0.05)

# Simulate data.
set.seed(1)
x <- rnorm(n)
y <- rnorm(n,x*b,s)

# Compute the posterior distribution of b.
xx    <- sum(x^2)
xy    <- sum(x*y)
bhat  <- xy/xx
shat  <- s/sqrt(xx)
sd    <- 1/sqrt(1/shat^2 + 1/s0^2)
mu    <- (sd/shat)^2*bhat
logBF <- log(shat/sqrt(shat^2 + s0^2)) + (mu/sd)^2/2
pp    <- w0*exp(logBF - max(logBF))
pp    <- pp/sum(pp)
