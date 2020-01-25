# Optimize the variational lower bound ("ELBO") for the mr-ash model
# by running a fixed number of co-ordinate ascent updates.
mr_ash <- function (X, y, b, se, s0, w0, numiter = 100) {

  # Center X and y.
  X <- scale(X,center = TRUE,scale = FALSE)
  y <- y - mean(y)
  
  # This variable is used to keep track of the algorithm's progress:
  # it stores the value of the objective (the variational lower bound,
  # or "ELBO") at each iteration.
  value <- rep(0,numiter)

  # Iterate the E and M steps.
  for (iter in 1:numiter) {

    # E STEP
    # ------
    # Run the co-ordinate ascent updates to update the variational
    # approximation to the posterior distribution of regression
    # coefficients.
    b <- mr_ash_update(X,y,b,se,s0,w)

    # M STEP
    # ------
    # Update the residual variance.
    # TO DO.
    
    # Update the mixture weights.
    # TO DO.
    
    # Record the algorithm's progress.
    value[iter] <- mr_ash_elbo(X,y,b,s,s0,w0)
  }

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = b,value = value))
}

# Do a single pass of the co-ordinate ascent updates.
mr_ash_update <- function (X, y, b, se, s0, w0) {

  # Get the number of predictors.
  p <- length(b)
  
  # Compute the expected residuals.
  r <- y - drop(X %*% b)
  
  # Repeat for each predictor.
  for (i in 1:p) {

    # Remove the jth effect from expected residuals.
    rj <- r + X[,j]*b[j]

    # Fit a Bayesian single-effect regression with a
    # mixture-of-normals prior to the expected residuals.
    out <- slr_mix(rj,X[,j],se,s0,w0)

    # Update the expected residuals.
    r <- rj - X[,j]*b[j]
  }

  # Return the updated posterior mean coefficients.
  return(b)
}

# Fit a single-effect regression model in which the regression
# coefficient is assigned a normal prior with mean zero and
# variance s0.
slr_ridge <- function (x, y, se, s0) {

  # Compute the least-squares estimate (bhat) and its variance (s).
  bhat <- dot(x,y)/dot(x,x)
  s    <- se/dot(x,x)

  # Compute the posterior mean (mu1) and variance (s1).
  s1  <- s0/(1 + s0/s)
  mu1 <- s1/s*bhat

  # Compute the log-Bayes factor.
  logbf <- (log(s/(s0 + s)) + (norm2(bhat)^2/s - norm2(bhat)^2/(s0 + s)))/2

  # Return the least-squares estimate (bhat), its variance (s), the
  # posterior mean and standard deviation (mu1, s1), and the
  # logarithm of the Bayes factor (logbf).
  return(list(bhat = bhat,s = s,mu1 = mu1,s1 = s1,logbf = logbf))
}
