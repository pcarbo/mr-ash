library(glmnet)
library(varbvs)
library(mr.ash.alpha)

# Analysis settings: number of simulations (ns), number of samples in
# training set (n), number of simulated variables (p), maximum number
# of variables with an effect on the continuous outcome (p1),
# proportion of variance in the outcome explained by the variables
# (pve).
ns  <- 100
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

# Predict the test set outcomes using the fitted models.
y.glmnet <- drop(predict(fit.glmnet,Xtest,s = "lambda.min"))

# Assess accuracy of the test predictions by compute the root-mean
# squared error (RMSE).
print(sqrt(mean((ytest - y.glmnet)^2)),digits = 3)

