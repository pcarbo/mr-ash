# TO DO: Explain here what this script does, and how to use it.
source("../code/misc.R")
source("../code/mr_ash.R")

# SCRIPT PARAMETERS
# ----------------
# Data simulation settings.
n  <- 500
p  <- 1000
sd <- c(0,    1,    2)
w  <- c(0.98, 0.01, 0.01)

# SIMULATE DATA
# -------------
set.seed(1)
X    <- matrix(rnorm(n*p),n,p)
k    <- sample(length(w),p,replace = TRUE,prob = w)
beta <- sd[k] * rnorm(p)
y    <- drop(X %*% beta + rnorm(n))

# FIT MR-ASH MODEL
# ----------------
b0  <- rep(0,p)
s   <- 1
s0  <- c(0.1,1,10)^2
w0  <- c(0.5,0.25,0.25)
fit <- mr_ash(X,y,s,s0,w0,b0,100)
plot(beta,fit$b,pch = 20,col = "black")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")
