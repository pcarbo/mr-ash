# Perform EM updates for the mr-ash model.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_ash <- function (X, y, se, s0, w0, b, numiter = 100) {

  # Center X and y.
  X <- scale(X,scale = FALSE)
  y <- y - mean(y)
  
  # These two variables are used to keep track of the algorithm's
  # progress: "elbo" stores the value of the objective (the
  # variational lower bound, or "ELBO") at each iteration; "maxd"
  # stores the largest different in the posterior mean coefficients
  # between two successive iterations.
  elbo <- rep(0,numiter)
  maxd <- rep(0,numiter)

  # Iterate the EM updates.
  for (i in 1:numiter) {

    # Save the current estimates of the posterior means.
    b0 <- b
    
    # E STEP
    # ------
    # Update the posterior means of the regression coefficients via
    # co-ordinate ascent.
    b <- mr_ash_update(X,y,b,se,s0,w)

    # M STEP
    # ------
    # Update the residual variance.
    # TO DO.
    
    # Update the mixture weights.
    # TO DO.
    
    # Record the algorithm's progress.
    # elbo[i] <- mr_ash_elbo(X,y,b,s,s0,w0)
    maxd[i] <- abs(max(b - b0))
  }

  # Return the updated posterior means of the regression coefficicents
  # ("b"), the value of the objective at each iteration ("elbo"), and
  # the maximum change at each iteration ("maxd").
  return(list(b    = b,
              elbo = elbo,
              maxd = maxd))
}

# Perform a single pass of the co-ordinate ascent updates for the
# mr-ash model.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_ash_update <- function (X, y, b, se, s0, w0) {

  # Get the number of predictors.
  p <- ncol(X)
  
  # Compute the expected residuals.
  r <- drop(y - X %*% b)
  
  # Repeat for each predictor.
  for (i in 1:p) {

    # Disregard the ith predictor in the expected residuals.
    r <- r + X[,i]*b[i]

    # Update the posterior distribution of the regression coefficients
    # for the ith predictor.
    out <- slr_mix(r,X[,i],se,s0,w0)

    # Update the expected residuals.
    r <- r - X[,i]*b[i]
  }

  # Return the updated posterior mean coefficients.
  return(b)
}

# Fit a univariate linear regression model in which the regression
# coefficient is assigned a normal prior with mean zero and
# variance s0.
bayes_lr_ridge <- function (x, y, se, s0) {

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
