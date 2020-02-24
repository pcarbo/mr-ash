# Fit the mr-ash model, in which the "naive" M-step update for the
# mixture weights is replaced by a mix-SQP update.
#
# This is a preliminary implementation and may not work in all cases.
mr_ash_with_mixsqp <- function (X, y, se, s0, w0, b, numiter = 10,
                                update.se = TRUE, maxiter.inner = 100,
                                tol.inner = 1e-8, verbose = TRUE) {

  # Get the number of predictors (p) and the number of mixture
  # components in the prior (k).
  p <- length(b)
  k <- length(w0)
  
  # Center X and y.
  X <- scale(X,scale = FALSE)
  y <- y - mean(y)

  # These two variables are used to keep track of the algorithm's
  # progress: "elbo" stores the value of the objective (the
  # variational lower bound, or "ELBO") at each iteration; "niter"
  # stores the number of inner-loop iterations performed at eacch
  # outer-loop iteration; and "maxd" stores the largest difference in
  # the mixture weights between two successive iterations.
  elbo  <- rep(0,numiter)
  maxd  <- rep(0,numiter)
  niter <- rep(0,numiter)
  
  # Iterate the mix-SQP updates.
  if (verbose) {
    cat("iter                elbo max|b-b'| max|w0-w0'| E-step M-step   ")
    cat("alpha\n")
  }
  for (i in 1:numiter) {

    # Save the current estimates of the mixture weights.
    w00 <- w0

    # Update the residual variance (se) and the posterior mean estimates
    # of the coefficients (b).
    out1 <- mr_ash(X,y,se,s0,w0,b,maxiter.inner,tol.inner,update.se,
                   update.w0 = FALSE,verbose = FALSE)
    se   <- out1$se
    b    <- out1$b

    # Update the mixture weights using the EM update.
    w0 <- out1$w0.em

    # Update the mixture weights using the mix-SQP update.
    out2 <- update_mixture_weights_with_mixsqp(X,y,se,s0,w0,b)
    w0   <- out2$w0

    # Record the algorithm's progress.
    elbo[i]  <- max(out1$elbo)
    niter[i] <- length(out1$elbo)
    maxd[i]  <- max(abs(w0 - w00))
    
    # Report progress, if requested.
    if (verbose)
      cat(sprintf("%4d %0.12e %0.3e %0.5e %6d %6d %0.1e\n",i,elbo[i],
                  tail(out1$maxd,n = 1),maxd[i],niter[i],out2$nmixsqp,
                  out2$a))
  }

  # Return the updated posterior means of the regression coefficicents
  # ("b"); the prior variances of the mixture components ("s0"); the
  # updated mixture weights ("w0"); the updated residual variance
  # ("se"); the value of the objective after running each outer-loop
  # iteration ("elbo"); the maximum change in the mixture weights at
  # each outer-loop iteration ("maxd"); and the number of inner-loop
  # (E-step) iterations performed for each outer-loop iteration
  # ("niter").
  return(list(b     = b,
              s0    = s0,
              w0    = w0,
              se    = se,
              elbo  = elbo,
              niter = niter,
              maxd  = maxd))
}

# TO DO: Explain here what this function does, and how to use it.
update_mixture_weights_with_mixsqp <- function (X, y, se, s0, w0, b,
                                                stepsize.reduce = 0.1,
                                                stepsize.min = 1e-8) {

  # Get the number of predictors (p) and the number of mixture
  # components in the prior (k).
  p <- length(b)
  k <- length(w0)

  # Compute the p x k matrix of log-likelihoods conditional on each
  # prior mixture component.
  L <- matrix(0,p,k)
  r <- drop(y - X %*% b)
  for (i in 1:p)
    for (j in 1:k)
      L[i,j] <- bayes_lr_ridge(X[,i],r + X[,i]*b[i],se,s0[j])$logbf
  out       <- mixsqp(L,log = TRUE,control = list(verbose = FALSE))
  w0.mixsqp <- out$x
  if (out$status != "converged to optimal solution")
    warning("mixsqp did not converge to optimal solution")

  # Compute the objective (ELBO) at the current iterate.
  f <- mr_ash_update(X,y,b,se,s0,w0)$elbo

  # Perform backtracking line search to identify a step size that
  # increases the ELBO.
  a   <- 1
  nls <- 0
  while (TRUE) {
    w0new <- a*w0.mixsqp + (1 - a)*w0
    fnew  <- mr_ash_update(X,y,b,se,s0,w0new)$elbo
    nls   <- nls + 1;

    # Check whether the new candidate increases the ELBO.
    if (fnew > f)
      break

    # If we cannot decrease the step size further, terminate the
    # backtracking line search, and set the step size to be the
    # minimum step size.
    else if (a*stepsize.reduce < stepsize.min) {

      # We need to terminate backtracking line search because we have
      # arrived at the smallest allowable step size.
      a     <- 0
      w0new <- w0
      break
    }

    # The new candidate does not increase the ELBO, so we need to try
    # again with a smaller step size.
    a = a * stepsize.reduce;
  }

  # Return the updated mixture weights (w0), the number of mix-SQP
  # iterations performed (nmixsqp), the number of backtracking line
  # search iterations (nls), and the step size determined by the
  # backtracking line search.
  return(list(w0      = w0new,
              nmixsqp = nrow(out$progress),
              nls     = nls,
              a       = a));
}
