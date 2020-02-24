# Perform EM updates for the mr-ash model.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_ash <- function (X, y, se, s0, w0, b, maxiter = 100, tol = 1e-8,
                    update.se = TRUE, update.w0 = TRUE, verbose = TRUE) {

  # Center X and y.
  X <- scale(X,scale = FALSE)
  y <- y - mean(y)

  # These two variables are used to keep track of the algorithm's
  # progress: "elbo" stores the value of the objective (the
  # variational lower bound, or "ELBO") at each iteration; "maxd"
  # stores the largest difference in the posterior mean coefficients
  # between two successive iterations.
  elbo <- rep(as.numeric(NA),maxiter)
  maxd <- rep(as.numeric(NA),maxiter)

  # Iterate the EM updates.
  if (verbose)
    cat("iter                elbo max|b-b'|\n")
  for (i in 1:maxiter) {

    # Save the current estimates of the posterior means.
    b0 <- b
    
    # E STEP
    # ------
    # Update the posterior means of the regression coefficients via
    # co-ordinate ascent.
    out   <- mr_ash_update(X,y,b,se,s0,w0)
    b     <- out$b
    w0.em <- out$w0.em
    
    # Record the algorithm's progress.
    elbo[i] <- out$elbo
    maxd[i] <- max(abs(b - b0))

    # M STEP
    # ------
    # Update the residual variance, if requested.
    if (update.se)
      se <- out$erss/n
    
    # Update the mixture weights, if requested.
    if (update.w0)
      w0 <- w0.em

    # Report progress, if requested.
    if (verbose)
      cat(sprintf("%4d %0.12e %0.3e\n",i,elbo[i],maxd[i]))

    # Stop if the largest change in the posterior means of the
    # regression coefficients is small.
    if (maxd[i] <= tol)
      break
  } 

  # Return the updated posterior means of the regression coefficicents
  # ("b"); the prior variances of the mixture components ("s0"); the
  # updated mixture weights ("w0"); the EM update for the mixture
  # weights ("w0.em"); residual variance ("se"); the value of the
  # objective after the (approximate) E-step at each iteration
  # ("elbo"); and the maximum change in the regression coefficients at
  # each iteration ("maxd").
  #
  # Note that w0 and w0.em will be the same whenever update.w0 = TRUE.
  return(list(b     = b,
              s0    = s0,
              w0    = w0,
              w0.em = w0.em,
              se    = se,
              elbo  = elbo[1:i],
              maxd  = maxd[1:i]))
}

# Perform a single pass of the co-ordinate ascent updates for the
# mr-ash model.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_ash_update <- function (X, y, b, se, s0, w0) {

  # Get the number of samples (n), the number of predictors (p), and
  # the number of mixture components (k).
  n <- nrow(X)
  p <- ncol(X)
  k <- length(w0)

  # This will be the M-step update for the weights in the
  # mixture-of-normals.
  w0.em <- rep(0,k)

  # Initialize two variables which will contain the expected residual
  # sum of squares (erss) and the variational lower bound (elbo) upon
  # completion of the loop. Also, store the term in the varaiationaal
  # lower bound that does not change as the variational approximation
  # is updated during the loop (elbo_const).
  elbo_const <- -n*log(2*pi*se)/2
  elbo       <- elbo_const
  erss       <- 0
  
  # Compute the expected residuals.
  r <- drop(y - X %*% b)

  # Repeat for each predictor.
  for (i in 1:p) {
    x <- X[,i]
      
    # Disregard the ith predictor in the expected residuals.
    r <- r + x*b[i]

    # Update the posterior distribution of the regression coefficients
    # for the ith predictor.
    out   <- bayes_lr_mix(x,r,se,s0,w0)
    b[i]  <- out$mu1
    w0.em <- w0.em + out$w1

    # Calculate the ith term in the "sum of variances" for the ELBO.
    v <- norm2(x)^2 * out$s1
    
    # Update the expected residuals.
    r <- r - x*b[i]
    
    # Compute the expected residual sum of of squares (erssi) for the
    # linear regression r = xi*bi + e, and the KL-divergence between
    # the posterior and prior distributions of the regression
    # coefficient for the ith predictor.
    erssi <- norm2(r)^2 + v
    d     <- elbo_const - (out$logbf + out$loglik0 + erssi/(2*se))
    elbo  <- elbo - d
    erss  <- erss + v
  }

  # Complete the computation of the expected residual sum of squares
  # (erss) and the variational lower bound (elbo).
  erss <- erss + norm2(r)^2
  elbo <- elbo - erss/(2*se)
  
  # Output the updated posterior mean coefficients ("b"), the M-step
  # update for the mixture weights ("w0.em"), the updated variational
  # lower bound (elbo), and the updated expected residual sum of
  # squares ("erss").
  return(list(b     = b,
              w0.em = w0.em/p,
              elbo  = elbo,
              erss  = erss))
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

  # Compute the "null" log-likelihood.
  loglik0 <- ldnorm(y,se)
  
  # Compute the log-Bayes factor.
  logbf <- ldnorm(bhat,s0 + s) - ldnorm(bhat,s)

  # Return the least-squares estimate (bhat) and its variance (s), the
  # posterior mean (mu1) and variance (S1), the log-likelihood for the
  # "null" model (when b = 0), and the log-Bayes factor (logbf).
  return(list(bhat    = bhat,
              s       = s,
              mu1     = mu1,
              s1      = s1,
              loglik0 = loglik0,
              logbf   = logbf))
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
  # posterior assignment probabilities (w1), the log-likelihood for
  # the "null" model (when b = 0), and the log-Bayes factor (logbf).
  return(list(mu1     = mu1,
              s1      = drop(s1),
              w1      = w1,
              loglik0 = out[[1]]$loglik0,
              logbf   = logbf))
}
