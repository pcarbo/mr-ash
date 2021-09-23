# TO DO: Explain here what this script is for, and how to use it.
library(varbvs)
library(mr.ash.alpha)

# This R code provides a very simple implementation of the mr-ash
# algorithm.
source("../code/misc.R")
source("../code/mr_ash.R")

# These are the data simulation settings.
n     <- 500
btrue <- c(1,-1,0,0)

# This specifies the variances for the mixture-of-normals prior on the
# regression coefficients.
s0 <- c(0.001,0.5,1)^2

# Simulate a data set.
set.seed(1)
p <- length(btrue)
X <- matrix(rnorm(n*p),n,p)
X <- scale(X,center = TRUE,scale = FALSE)
y <- drop(X %*% btrue + rnorm(n))
y <- y - mean(y)

# These are the initial estimates of residual variance (s), mixture
# weights (w0) and posterior mean estimates of the regression
# coefficients (b).
s  <- 1
w0 <- rep(1/3,3)
b  <- rep(0,p)

# Fit the model using coordinate ascent updates.
fit1 <- mr_ash(X,y,s,s0,w0,b,method = "cd",maxiter = 20)

# Fit the model using Nelder-Mead.
fit2 <- mr_ash(X,y,s,s0,w0,b,method = "nm",maxiter = 20)
           
# Check my calculations against mr.ash.alpha and varbvsmix.
fit3 <- with(fit1,
             mr.ash(X,y,sa2 = c(0,s0[-1])/se,sigma2 = se,pi = w0,
                    update.pi = FALSE,update.sigma2 = FALSE,
                    max.iter = 20))
fit4 <- with(fit1,
             varbvsmix(X,NULL,y,c(0,s0[-1])/se,se,w0,
                  update.sigma = FALSE,update.w = FALSE,
                  maxiter = 20,tol = 1e-8,verbose = FALSE))
b1 <- with(fit4,rowSums(alpha*mu))
cat(sprintf("mr_ash(method=cd): %0.12f\n",max(fit1$elbo)))
cat(sprintf("mr_ash(method=nm): %0.12f\n",max(fit2$elbo)))
cat(sprintf("mr.ash.alpha:      %0.12f\n",max(-fit3$varobj)))
cat(sprintf("varbvsmix:         %0.12f\n",max(fit4$logZ + log(n)/2)))
