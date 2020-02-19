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
plot(max(fit$elbo) - fit$elbo + 1e-8,type = "l",col = "dodgerblue",
     lwd = 2,log = "y",xlab = "iteration",ylab = "distance to best ELBO")

# TESTING
# -------
library(mixsqp)
numiter <- 20
k    <- 3
se   <- s
X    <- scale(X,scale = FALSE)
elbo <- rep(0,numiter)
for (iter in 1:numiter) {

  # Update residual variance and posterior on coefficients.
  fit <- mr_ash(X,y,se,s0,w0,b,100,tol = 1e-6,update.w0 = FALSE)
  se  <- fit$se
  b   <- fit$b
  elbo[iter] <- max(fit$elbo)
  
  # Update mixture weights using mixsqp.
  L <- matrix(0,p,k)
  r <- drop(y - X %*% b)
  for (i in 1:p) {
    x <- X[,i]
    r <- r + x*b[i]
    for (j in 1:k)
      L[i,j] <- bayes_lr_ridge(x,r,se,s0[j])$logbf
    r <- r - x*b[i]
  }
  out <- mixsqp(L,log = TRUE,control = list(verbose = FALSE))
  w0  <- out$x
}
