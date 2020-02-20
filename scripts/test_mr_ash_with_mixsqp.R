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

scenario <- "independent"

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
fit1 <- mr_ash(X,y,s,s0,w0,b,100)

# Fit the model a second time using the "accelerated" mix-SQP updates
# for the mixture weights.
fit2 <- mr_ash_with_mixsqp(X,y,s,s0,s0,b,10)

stop()

# REVIEW FITS
# -----------
# Compare the posterior mean estimates against the values used to
# simulate the data.
plot(beta,fit$b,pch = 20,col = "black")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

# Plot the improvement in the solution over time.
plot(max(fit$elbo) - fit$elbo + 1e-8,type = "l",col = "dodgerblue",
     lwd = 2,log = "y",xlab = "iteration",ylab = "distance to best ELBO")
