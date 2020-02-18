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
    out <- mr_ash_update(X,y,b,se,s0,w0)
    b   <- out$b
    
    # Record the algorithm's progress.
    elbo[i] <- out$elbo
    maxd[i] <- abs(max(b - b0))
    
    # M STEP
    # ------
    # Update the residual variance.
    # TO DO.
    
    # Update the mixture weights.
    w0 <- out$w0.em
  } 

  # Return the updated posterior means of the regression coefficicents
  # ("b"), the updated mixture weights ("w0"), the value of the
  # objective after the (approximate) E-step at each iteration
  # ("elbo"), and the maximum change at each iteration ("maxd").
  return(list(b    = b,
              w0   = w0,
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

  # Get the number of samples (n), the number of predictors (p), and
  # the number of mixture components (k).
  n <- nrow(X)
  p <- ncol(X)
  k <- length(w0)

  # FOR TESTING ONLY
  mu1_mix <- matrix(0,p,k)
  s1_mix <- matrix(0,p,k)
  w1_mix <- matrix(0,p,k)
  
  # This will be the M-step update for the weights in the
  # mixture-of-normals.
  w0.em <- rep(0,k)

  # This is used to store the store the sum-of-variances term in the
  # expression for the ELBO.
  v <- 0
  d <- 0
  
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
    v     <- v + norm2(x)^2 * out$s1
    
    # FOR TESTING ONLY
    mu1_mix[i,] <- out$mu1_mix
    s1_mix[i,]  <- out$s1_mix
    w1_mix[i,]  <- out$w1

    f <- 0
    s1 <- out$s1_mix
    mu1 <- out$mu1_mix
    w1 <- out$w1
    for (j in 1:k) {
      f <- f + w1[j]*log(w0[j]) - w1[j]*log(w1[j])
      f <- f + w1[j]/2 + w1[j]*log(s1[j]/s0[j])/2 -
               (w1[j]*(s1[j] + mu1[j]^2))/(2*s0[j])
    }
    d <- d - f
    
    # Update the expected residuals.
    r <- r - x*b[i]
  }

  # FOR TESTING ONLY
  # elbo  <- compute_elbo(X,y,se,s0,w0,w1_mix,mu1_mix,s1_mix)
  erss <- norm2(r)^2
  elbo <- -n/2*log(2*pi*se) - (erss + v)/(2*se) - d
  
  # Output the updated posterior mean coefficients (b) and the M-step
  # update for the mixture weights (w0.em).
  return(list(b     = b,
              w0.em = w0.em/p,
              elbo  = elbo))
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
              logbf = logbf,

              # FOR TESTING ONLY
              mu1_mix = sapply(out,"[[","mu1"),
              s1_mix  = sapply(out,"[[","s1")))
}

# FOR TESTING ONLY.
compute_elbo <- function (X, y, se, s0, w0, w1, mu1, s1) {
  e <- 1e-30
  
  # Get the number of samples (n), the number of variables (p), and
  # the number of mixture components (k).
  n <- nrow(X)
  p <- ncol(X)
  k <- length(w0)

  # Compute the ELBO.
  d <- diag(crossprod(X))
  b <- rowSums(w1 * mu1)
  f <- -n*log(2*pi*se)/2 -
        norm2(y - X %*% b)^2/(2*se) - 
        dot(d,betavarmix(w1,mu1,s1))/(2*se)
  for (i in 1:k)
    f <- f + sum(w1[,i]*log(w0[i] + e)) - sum(w1[,i]*log(w1[,i] + e))
  for (i in 1:k)
    f <- f + sum(w1[,i])/2 + dot(w1[,i],log(s1[,i]/s0[i]))/2 -
             dot(w1[,i],s1[,i] + mu1[,i]^2)/(2*s0[i])
  return(f)
}

# FOR TESTING ONLY.
betavarmix <- function (p, mu, s)
  rowSums(p*(s + mu^2)) - rowSums(p*mu)^2
