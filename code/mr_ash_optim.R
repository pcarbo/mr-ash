# TO DO: Explain here what this function does, and how to use it.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_ash_optim <- function (X, y, se, s0, w0, b, method = c("Nelder-Mead"),
                          maxiter = 100, tol = 1e-8, update.se = TRUE,
                          update.w0 = TRUE, verbose = TRUE) {

  # Center X and y.
  X <- scale(X,scale = FALSE)
  y <- y - mean(y)

  # Compute theta.
  xx    <- diag(crossprod(X))
  theta <- drop(y %*% X)/xx
  
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
    # Update the posterior means of the regression coefficients using
    # optim.
    #
    # *** TO DO ***
    #
    optim()
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
