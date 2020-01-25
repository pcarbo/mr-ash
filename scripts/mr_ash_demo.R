# TO DO: Explain here what this script does, and how to use it.
source("../code/misc.R")
source("../code/mr_ash.R")

# SCRIPT PARAMETERS
# ----------------
# Data simulation settings.
n    <- 500
p    <- 1000
sd   <- c(0,    1,    2)
w    <- c(0.98, 0.01, 0.01)

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
s0  <- c(0.1,1,10)^2
s   <- 1
w   <- c(0.5,0.25,0.25)
fit <- mr_ash(X,y,b0,s,s0,w,100)
