library(MASS)
library(mixsqp)
source("../code/misc.R")
source("../code/mr_ash.R")

# Script settings.
n  <- 100
p  <- 400
sd <- c(0,   1,    2)
w  <- c(0.9, 0.05, 0.05)
s  <- 0.1
s0 <- 10^seq(-4,0,length.out = 12)

# Simulate data.
set.seed(2)
S       <- matrix(s,p,p)
diag(S) <- 1
X       <- mvrnorm(n,rep(0,p),S)
k       <- sample(length(w),p,replace = TRUE,prob = w)
beta    <- sd[k] * rnorm(p)
y       <- drop(X %*% beta + rnorm(n))
k  <- length(s0)
se <- 1
w0 <- rep(1/k,k)
b  <- rep(0,p)

# Fit model.
fit1 <- mr_ash(X,y,se,s0,w0,b,maxiter = 200,verbose = TRUE)
fit2 <- mr_ash_with_mixsqp(X,y,se = fit1$se,s0,w0,b,numiter = 10)

# Plot improvement in solution over time.
elbo.best <- max(c(fit1$elbo,fit2$elbo))
pdat      <- rbind(data.frame(update = "em",
                              iter   = 1:length(fit1$elbo),
                              elbo   = fit1$elbo),
                   data.frame(update = "mixsqp",
                              iter   = cumsum(fit2$niter),
                              elbo   = fit2$elbo))
pdat$elbo <- elbo.best - pdat$elbo + 1e-4
ggplot(pdat,aes(x = iter,y = elbo,color = update)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  scale_color_manual(values = c("royalblue","darkorange")) +
  labs(y = "distance to \"best\" elbo") +
  theme_cowplot()
