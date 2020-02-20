# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(cowplot)
library(MASS)
library(mixsqp)
source("../code/misc.R")
source("../code/mr_ash.R")

# SCRIPT PARAMETERS
# ----------------
# Data simulation settings.
n  <- 500
p  <- 100
sd <- c(0,   1,    2)
w  <- c(0.9, 0.05, 0.05)
s  <- 0.5

# Variances for the mixture-of-normals prior on the regression
# coefficients.
s0 <- 10^seq(-4,0,length.out = 12)

# SIMULATE DATA
# -------------
# The predictors are drawn from the multivariate normal with zero mean
# and covariance matrix S, in which all diagonal entries are 1, and
# all off-diagonal entries are s. Setting s = 0.5 reproduces the
# simulation of the predictors used in Example 3 of Zou & Hastie
# (2005).
cat("Simulating data.\n")
set.seed(2)
S       <- matrix(s,p,p)
diag(S) <- 1
X       <- mvrnorm(n,rep(0,p),S)
k       <- sample(length(w),p,replace = TRUE,prob = w)
beta    <- sd[k] * rnorm(p)
y       <- drop(X %*% beta + rnorm(n))

# FIT MR-ASH MODEL
# ----------------
# These are the initial estimates of residual variance (s), mixture
# weights (w0), and posterior mean estimates of the regression
# coefficients (b).
k   <- length(s0)
s   <- 1
w0  <- rep(1/k,k)
b   <- rep(0,p)

# Fit the model by running 100 EM updates.
cat("Running EM updates.\n")
fit1 <- mr_ash(X,y,s,s0,w0,b,maxiter = 500,verbose = FALSE)

# Fit the model a second time using the "accelerated" mix-SQP updates
# for the mixture weights.
cat("Running 10 mix-SQP updates.\n")
fit2 <- mr_ash_with_mixsqp(X,y,s,s0,w0,b,numiter = 10,tol.inner = 1e-4)

# REVIEW FITS
# -----------
# Compare the posterior mean estimates against the values used to
# simulate the data.
cat("Summarizing results.\n")
p1 <- quickplot(beta,fit1$b) +
  geom_abline(intercept = 0,slope = 1,col = "skyblue",lty = "dotted") +
  theme_cowplot()

# Compare the posterior mean estimates from the two fits.
p2 <- quickplot(fit1$b,fit2$b) +
  geom_abline(intercept = 0,slope = 1,col = "skyblue",lty = "dotted") +
  theme_cowplot()

# Plot the improvement in the solution over time.
elbo.best <- max(c(fit1$elbo,fit2$elbo))
pdat      <- rbind(data.frame(update = "em",
                              iter   = 1:length(fit1$elbo),
                              elbo   = fit1$elbo),
                   data.frame(update = "mixsqp",
                              iter   = cumsum(fit2$niter),
                              elbo   = fit2$elbo))
pdat$elbo <- elbo.best - pdat$elbo + 1e-8
p3 <- ggplot(pdat,aes(x = iter,y = elbo,color = update)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  scale_color_manual(values = c("darkblue","darkorange")) +
  theme_cowplot()
