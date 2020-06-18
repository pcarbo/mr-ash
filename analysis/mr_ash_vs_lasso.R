library(mixsqp)
library(glmnet)
library(mr.ash.alpha)
library(ggplot2)
library(cowplot)
source("../code/functions_for_mr_ash_vs_lasso.R")

set.seed(1)

# SCRIPT SETTINGS
# ---------------
ns  <- 10   # Number of simulations.
n   <- 500  # Number of samples in training set.
p   <- 1000 # Number of candidate predictors.
p1  <- 200  # Number of predictors with nonzero effects.
pve <- 0.95 # Proportion of variance explained.

# These two vectors will store the results from each simulation.
rmse1 <- rep(0,ns)
rmse2 <- rep(0,ns)
rmse3 <- rep(0,ns)
rmse4 <- rep(0,ns)

# Repeat for each simulation.
for (i in 1:ns) {
  cat(i,"")
  
  # SIMULATE DATA
  # -------------
  X    <- matrix(rnorm(2*n*p),2*n,p)
  X    <- scale(X,center = TRUE,scale = TRUE)
  b    <- rep(0,p)
  j    <- sample(p,p1)
  b[j] <- rnorm(p1)
  y    <- drop(X %*% b)
  se   <- sqrt((1 - pve)/pve)*sd(y)
  y    <- y + rnorm(n,sd = se)
  y    <- y - mean(y)

  # Split the data 50-50 into a training set and a test set.
  test  <- 1:n
  Xtest <- X[test,]
  ytest <- y[test]
  X     <- X[-test,]
  y     <- y[-test]

  # FIT LINEAR REGRESSIONS
  # ----------------------
  # Run Lasso model, in which the penalty strength is estimated by
  # 10-fold cross-validation.
  fit1 <- cv.glmnet(X,y,standardize = FALSE)

  # Run mr.ash.
  fit3  <- mr.ash(X,y,standardize = FALSE)

  # Run mr.ash a second time, in which the residual variance parameter
  # and the regression coefficients are initialized from the Lasso estimates.
  b    <- coef(fit1,s = "lambda.min")[-1]
  yest <- drop(predict(fit1,X,s = "lambda.min"))
  s    <- var(y - yest)
  fit4 <- mr.ash(X,y,beta.init = b,sigma2 = s)

  # PREDICT TEST OUTCOMES
  # ---------------------
  # Compute fitted values in test set.
  y1 <- drop(predict(fit1,Xtest,s = "lambda.min"))
  y2 <- drop(predict(fit1,Xtest,s = "lambda.1se"))
  y3 <- predict(fit3,Xtest)
  y4 <- predict(fit4,Xtest)

  # Assess the accuracy of the predictions.
  rmse1[i] <- sqrt(mean((ytest - y1)^2))
  rmse2[i] <- sqrt(mean((ytest - y2)^2))
  rmse3[i] <- sqrt(mean((ytest - y3)^2))
  rmse4[i] <- sqrt(mean((ytest - y4)^2))
  print(rmse1[i])
  print(rmse2[i])
  print(rmse3[i])
  print(rmse4[i])
  stop()
}
cat("\n")

# SUMMARIZE RESULTS
# -----------------
# Create a scatterplot comparing the accuracy of the mr.ash estimates
# against the glmnet estimates.
p1 <- ggplot(data.frame(glmnet = rmse1,mr.ash = rmse2),
       aes(x = glmnet,y = mr.ash)) +
  geom_point(size = 2,shape = 21,color = "white",fill = "black") +
  geom_abline(slope = 1,intercept = 0,color = "skyblue",linetype = "dotted") +
  xlim(16,22) +
  ylim(16,22) +
  theme_cowplot()

# Create a scatterplot comparing the accuracy of the mr.ash estimates,
# when fixing the residual variance parameter, against the glmnet
# estimates.
p2 <- ggplot(data.frame(glmnet = rmse1,mr.ash = rmse3),
       aes(x = glmnet,y = mr.ash)) +
  geom_point(size = 2,shape = 21,color = "white",fill = "black") +
  geom_abline(slope = 1,intercept = 0,color = "skyblue",linetype = "dotted") +
  xlim(16,25) +
  ylim(16,25) +
  labs(y = "mr.ash (true \u03c3)") +
  theme_cowplot(12)
