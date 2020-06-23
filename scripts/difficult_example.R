library(glmnet)
library(varbvs)
library(mr.ash.alpha)

# Analysis settings: number of simulations (ns), number of samples in
# training set (n), number of simulated variables (p), maximum number
# of variables with an effect on the continuous outcome (p1),
# proportion of variance in the outcome explained by the variables
# (pve).
n   <- 500
p   <- 1000
p1  <- 467
pve <- 0.95

# Simulate data.
set.seed(15)

# Simulate a 2n x p design matrix; the first n rows will be training
# data, and the last n rows will be test data.
X <- matrix(rnorm(2*n*p),2*n,p)
X <- scale(X,center = TRUE,scale = TRUE)

# Simulate the p regression coefficients; only p1 < p of the
# coefficients are nonzero.
b    <- rep(0,p)
j    <- sample(p,p1)
b[j] <- rnorm(p1)

# Simulate the responses.
y  <- drop(X %*% b)
se <- sqrt((1 - pve)/pve) * sd(y)
y  <- y + rnorm(n,sd = se)

# Split the data 50-50 into a training set and a test set.
test  <- 1:n
Xtest <- X[test,]
ytest <- y[test]
X     <- X[-test,]
y     <- y[-test]

# Fit the Elastic Net model, in which the penalty strength parameter
# ("lambda") is chosen via 10-fold cross-validation.
fit.glmnet <- cv.glmnet(X,y,alpha = 0.95,standardize = FALSE)

# Fit the mr.ash model. 
fit.mrash <- mr.ash(X,y,standardize = FALSE)

# Fit the mr.ash model a second time, but provide an (unfair)
# advantage to mr.ash by providing it with the prior used to simulate the
# data.
w1 <- p1/p
s  <- se^2
fit.optg <- mr.ash(X,y,beta.init = coef(fit.glmnet)[-1],
                   update.pi = FALSE,update.sigma2 = FALSE,
                   sigma2 = s,sa2 = c(0,1/s),pi = c(1 - w1,w1))

# Now we provide mr.ash with a slightly weaker (unfair) advantage,
# where we initialize it to the prior used to simulate the data, but
# allow it to fit the prior.
fit.initoptg <- mr.ash(X,y,beta.init = coef(fit.glmnet)[-1],
                       update.pi = TRUE,update.sigma2 = FALSE,
                       sigma2 = s,sa2 = c(0,1/s),pi = c(1 - w1,w1))

fit.varbvs <- varbvs(X,NULL,y,update.sigma = FALSE,sa = 1/s,
                     logodds = seq(-3,-1,length.out = 50),
                     verbose = FALSE)

# elbo <- rep(0,50)
# w1   <- 10^seq(-3,-1,length.out = 50)
# for (i in 1:50) {
#   fit <- mr.ash(X,y,beta.init = coef(fit.glmnet)[-1],update.pi = FALSE,
#                 update.sigma2 = FALSE,sigma2 = s,sa2 = c(0,1/s),
#                 pi = c(1 - w1[i],w1[i]))
#   elbo[i] <- -min(fit$varobj)
# }
# elbo <- elbo - max(elbo)

# Predict the test set outcomes using the fitted models.
y.glmnet   <- drop(predict(fit.glmnet,Xtest,s = "lambda.min"))
y.mrash    <- predict(fit.mrash,Xtest)
y.optg     <- predict(fit.optg,Xtest)
y.initoptg <- predict(fit.initoptg,Xtest)
y.varbvs   <- predict(fit.varbvs,Xtest)

# Assess accuracy of the test predictions by compute the root-mean
# squared error (RMSE).
print(sqrt(mean((ytest - y.glmnet)^2)),digits = 3)
print(sqrt(mean((ytest - y.mrash)^2)),digits = 3)
print(sqrt(mean((ytest - y.optg)^2)),digits = 3)
print(sqrt(mean((ytest - y.initoptg)^2)),digits = 3)
print(sqrt(mean((ytest - y.varbvs)^2)),digits = 3)

sigmoid10 <- function (x) 1/(1 + 10^(-x))
logodds <- fit.varbvs$logodds
logw <- fit.varbvs$logw
logw <- logw - max(logw)
plot(sigmoid10(logodds),logw,type="l",lwd = 2,col = "darkblue",log = "x")
points(sigmoid10(logodds),logw,pch = 20,col = "darkblue")

