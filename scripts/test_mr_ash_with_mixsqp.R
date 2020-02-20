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
#
# Suggestion: Try setting "s" to zero (all predictors are independent)
# and 0.4 (correlation between all predictors is 0.4).
#
n  <- 500
p  <- 100
sd <- c(0,   1,    2)
w  <- c(0.9, 0.05, 0.05)
s  <- 0.4

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
se  <- 1
w0  <- rep(1/k,k)
b   <- rep(0,p)

# Fit the model by running 500 EM updates.
cat("Running EM updates.\n")
fit1 <- mr_ash(X,y,se,s0,w0,b,maxiter = 500,verbose = FALSE)

# Fit the model a second time using the "accelerated" mix-SQP updates
# for the mixture weights.
cat("Running 10 mix-SQP updates.\n")
fit2 <- mr_ash_with_mixsqp(X,y,se,s0,w0,b,numiter = 10,tol.inner = 1e-8)

# REVIEW FITS
# -----------
# Compare the posterior mean estimates against the values used to
# simulate the data.
cat("Summarizing results.\n")
p1 <- ggplot(data.frame(true = beta,mixsqp = fit2$b),
             aes(x = true,y = mixsqp)) +
  geom_point(color = "darkblue") +
  geom_abline(intercept = 0,slope = 1,col = "magenta",lty = "dotted") +
  labs(title = "coefs") +
  theme_cowplot(10)

# Compare the posterior mean estimates from the two fits.
p2 <- ggplot(data.frame(em = fit1$b,mixsqp = fit2$b),
             aes(x = em,y = mixsqp)) +
  geom_point(color = "darkblue") +
  geom_abline(intercept = 0,slope = 1,col = "magenta",lty = "dotted") +
  labs(x = "em",y = "mixsqp",title = "coefs") +
  theme_cowplot(10)

# Plot the improvement in the solution over time. The EM updates are
# shown in blue, and the mix-SQP updates are shown in orange.
elbo.best <- max(c(fit1$elbo,fit2$elbo))
pdat      <- rbind(data.frame(update = "em",
                              iter   = 1:length(fit1$elbo),
                              elbo   = fit1$elbo),
                   data.frame(update = "mixsqp",
                              iter   = cumsum(fit2$niter),
                              elbo   = fit2$elbo))
pdat$elbo <- elbo.best - pdat$elbo + 0.01
p3 <- ggplot(pdat,aes(x = iter,y = elbo,color = update)) +
  geom_line(show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  scale_y_log10() +
  scale_color_manual(values = c("royalblue","darkorange")) +
  labs(x = "iteration",y = "distance to \"best\" ELBO") +
  theme_cowplot(10)

# Show all the plots.
print(plot_grid(p1,p2,p3,ncol = 3,rel_widths = c(2,2,2.5)))
