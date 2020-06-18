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
p1  <- 500  # Number of predictors with nonzero effects.
pve <- 0.95 # Proportion of variance explained.

# These two vectors will store the results from each simulation.
rmse1 <- rep(0,ns)
rmse2 <- rep(0,ns)

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

  # Split the data into training and test sets.
  test  <- 1:500
  Xtest <- X[test,]
  ytest <- y[test]
  X     <- X[-test,]
  y     <- y[-test]

  # FIT LASSO
  # ---------
  # Fit the Lasso model.
  fit1 <- cv.glmnet(X,y,standardize = FALSE)

  # FIT MR.ASH
  # ----------
  # Fit the mr.ash model.
  fit2 <- mr.ash(X,y)

  # PREDICT TEST OUTCOMES
  # ---------------------
  # Compute fitted values in test set.
  y1 <- predict(fit1,Xtest,s = "lambda.min")
  y2 <- predict(fit2,Xtest)

  # Assess the accuracy of the predictions.
  rmse1[i] = sqrt(mean((ytest - y1)^2))
  rmse2[i] = sqrt(mean((ytest - y2)^2))
}
cat("\n")

# SUMMARIZE RESULTS
# -----------------
# Create a scatterplot comparing the accuracy of the mr.ash estimates
# against the glmnet estimates.
ggplot(data.frame(glmnet = rmse1,mr.ash = rmse2),
       aes(x = glmnet,y = mr.ash)) +
  geom_point(size = 2,shape = 21,color = "white",fill = "black") +
  geom_abline(slope = 1,intercept = 0,color = "skyblue",linetype = "dotted") +
  xlim(range(c(rmse1,rmse2))) +
  ylim(range(c(rmse1,rmse2))) +
  theme_cowplot()
