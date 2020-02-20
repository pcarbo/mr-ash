# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(cowplot)
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

# REVIEW FITS
# -----------
# Compare the posterior mean estimates against the values used to
# simulate the data.
quickplot(beta,fit1$b) +
  geom_abline(intercept = 0,slope = 1,col = "skyblue",lty = "dotted") +
  theme_cowplot()

# Compare the posterior mean estimates from the two fits.
quickplot(fit1$b,fit2$b) +
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
ggplot(pdat,aes(x = iter,y = elbo,color = update)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  scale_color_manual(values = c("darkblue","darkorange")) +
  theme_cowplot()
