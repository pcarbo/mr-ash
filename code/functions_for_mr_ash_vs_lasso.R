# TO DO: Explain here what this function does, and how to use it.
simulate_data <- function (n, p, p1, pve) {

  # Simulate a 2n x p design matrix; the first n rows will be training
  # data, and the last n rows will be test data.
  X <- matrix(rnorm(2*n*p),2*n,p)
  X <- scale(X,center = TRUE,scale = TRUE)

  # Simulate the p regression coefficients; only p1 < p of the
  # coefficients are nonzero.
  b    <- rep(0,p)
  j    <- sample(p,p1)
  b[j] <- rnorm(p1)

  # Simulate the responses.
  y  <- drop(X %*% b)
  se <- sqrt((1 - pve)/pve) * sd(y)
  y  <- y + rnorm(n,sd = se)
  y  <- y - mean(y)

  # Split the data 50-50 into a training set and a test set.
  test  <- 1:n
  Xtest <- X[test,]
  ytest <- y[test]
  X     <- X[-test,]
  y     <- y[-test]

  # Output ...
  return(list(X     = X,
              y     = y,
              Xtest = Xtest,
              ytest = ytest,
              b     = b,
              se    = se))
}

