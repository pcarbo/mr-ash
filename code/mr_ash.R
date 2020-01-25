# Optimize the variational lower bound ("ELBO") for the mr-ash model
# by running a fixed number of co-ordinate ascent updates.
mr_ash <- function (X, y, b0, s, s0, w, numiter = 100) {

  # Get the initial estimates of the posterior mean regression
  # coefficients.
  b <- b0

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
    b <- mr_ash_update(X,y,b,s,s0,w)

    # M STEP
    # ------
    # Update the residual variance.
    # TO DO.
    
    # Update the mixture weights.
    # TO DO.
    
    # Record the algorithm's progress.
    value[iter] <- mr_ash_elbo(X,y,b,s,s0,w)
  }

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = b,value = value))
}

# Do a single pass of the co-ordinate ascent updates.
mr_ash_update <- function (X, y, b, s, s0, w) {

  # Get the number of variables (p) and the number of mixture
  # components (k).
  p <- length(b)
  k <- length(w)
  
  # Compute the expected residuals.
  r <- y - drop(X %*% b)
  
  # Repeat for each variable.
  for (i in 1:p) {

    # Remove the jth effect from expected residuals.
    rj <- r + X[,j]*b[j]

    # Fit a Bayesian single-effect regression.
    # TO DO.

    # Update the expected residuals.
    r <- rj - X[,j]*b[j]
  }
  
  return(b)
}
