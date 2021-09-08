# TO DO: Explain here what this script does, and how to use it.

# Script settings.
n  <- 75
s  <- 3
b  <- -1
s0 <- c(0.01,0.1,1)
w0 <- c(0.89,0.05,0.05)

# Simulate data.
set.seed(1)
x <- rnorm(n)
y <- rnorm(n,x*b,s)

# Compute the posterior distribution of b.
xx   <- sum(x^2)
xy   <- sum(x*y)
bhat <- xy/xx
shat <- s/sqrt(xx)
s1   <- 1/sqrt(1/shat^2 + 1/s0^2)
mu1  <- (s1/shat)^2*bhat
BF   <- shat/sqrt(shat^2 + s0^2) * exp((mu1/s1)^2/2)
p1   <- w0*BF
p1   <- p1/sum(p1)

# Compute the KL-divergence term.
kl1 <- -(sum(p1*(1 + log((s1/s0)^2) - (mu1^2 + s1^2)/s0^2))/2 +
         sum(p1*log(w0/p1)))

# Compute the KL-divergence term using the alternative expression.
b1   <- sum(p1*mu1)
v1   <- sum(p1*(mu1^2 + s1^2)) - b1^2
erss <- sum((y - x*b1)^2) + xx*v1
kl2  <- -log(sum(w0*BF)) - (b1^2 + v1 - 2*bhat*b1)/(2*shat^2)

# Compare the two calculations.
print(kl1 - kl2)
