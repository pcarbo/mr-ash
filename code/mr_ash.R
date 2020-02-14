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
    b <- mr_ash_update(X,y,b,se,s0,w0)

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
    x <- X[,i]
      
    # Disregard the ith predictor in the expected residuals.
    r <- r + x*b[i]

    # Update the posterior distribution of the regression coefficients
    # for the ith predictor.
    out  <- bayes_lr_mix(x,r,se,s0,w0)
    b[i] <- out$mu1

    # Update the expected residuals.
    r <- r - x*b[i]
  }

  # Output the updated posterior mean coefficients.
  return(b)
}

# Fit a univariate linear regression model in which the regression
# coefficient is assigned a normal prior with mean zero and
# variance s0.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
bayes_lr_ridge <- function (x, y, se, s0) {

  # Compute the least-squares estimate of the coefficients (bhat) and
  # its standard error (s).
  xx   <- norm2(x)^2
  bhat <- dot(x,y)/xx
  s    <- se/xx
  
  # Compute the posterior mean (mu1) and variance (s1) assuming a
  # normal prior wiith zero mean and variance s0.
  s1  <- s0/(1 + s0/s)
  mu1 <- s1/s*bhat

  # Compute the log-Bayes factor.
  logbf <- ldnorm(bhat,s0 + s) - ldnorm(bhat,s)

  # Return the least-squares estimate (bhat) and its variance (s), the
  # posterior mean (mu1) and variance (S1), and the log-Bayes factor
  # (logbf).
  return(list(bhat  = bhat,
              s     = s,
              mu1   = mu1,
              s1    = s1,
              logbf = logbf))
}

# Fit a univariate linear regression model in which the regression
# coefficient is assigned a mixture-of-normals prior.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
bayes_lr_mix <- function (x, y, se, s0, w0) {

  # Get the number of mixture components.
  k <- length(w0)

  # Compute the quantities separately for each mixture component.
  out <- vector("list",k)
  for (i in 1:k)
    out[[i]] <- bayes_lr_ridge(x,y,se,s0[i])

  # Compute the posterior assignment probabilities for the latent
  # indicator variable.
  logbf <- sapply(out,"[[","logbf")
  w1    <- softmax(logbf + log(w0))

  # Compute the log-Bayes factor as a linear combination of the
  # individual Bayes factors for each mixture component.
  z     <- log(w0) + logbf
  u     <- max(z)
  logbf <- u + log(sum(exp(z - u)))
  
  # Compute the posterior mean (mu1) and variance (s1) of the
  # regression coefficients.
  s1  <- 0
  mu1 <- 0
  for (i in 1:k) {
    w   <- w1[i]
    mu  <- out[[i]]$mu1
    s   <- out[[i]]$s1
    mu1 <- mu1 + w*mu
    s1  <- s1 + w*(s + mu^2)
  }
  s1 <- s1 - mu1^2
  
  # Return the the posterior mean (mu1) and svariance (s1), the
  # posterior assignment probabilities (w1), and the log-Bayes factor
  # (logbf).
  return(list(mu1   = mu1,
              s1    = drop(s1),
              w1    = w1,
              logbf = logbf))
}

