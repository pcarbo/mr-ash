# An illustration of the mr_ash algorithm applied to a small,
# simulated data set.
source("../code/misc.R")
source("../code/mr_ash.R")

# SCRIPT PARAMETERS
# ----------------
# Data simulation settings.
n  <- 500
p  <- 1000
sd <- c(0,    1,    2)
w  <- c(0.98, 0.01, 0.01)

# Variances for the mixture-of-normals prior on the regression
# coefficients.
s0 <- c(0.01,0.5,1)^2

# SIMULATE DATA
# -------------
set.seed(1)
X    <- matrix(rnorm(n*p),n,p)
k    <- sample(length(w),p,replace = TRUE,prob = w)
beta <- sd[k] * rnorm(p)
y    <- drop(X %*% beta + rnorm(n))

# FIT MR-ASH MODEL
# ----------------
# These are the initial estimates of residual variance (s), mixture
# weights (w0) and posterior mean estimates of the regression
# coefficients (b).
b   <- rep(0,p)
s   <- 1
w0  <- c(0.5,0.25,0.25)

# Fit the model by running 100 EM updates.
fit <- mr_ash(X,y,s,s0,w0,b,100)

# Compare the posterior mean estimates against the values used to
# simulate the data.
plot(beta,fit$b,pch = 20,col = "black")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

# Plot the improvement in the solution over time.
# TO DO.
