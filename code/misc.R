# Returns the dot product of vectors x and y.
dot <- function (x,y)
  sum(x*y)

# Returns the quadratic norm (2-norm) of vector x.
norm2 <- function (x)
  sqrt(dot(x,x))

# Compute the softmax of vector x in a more numerically prudent manner
# that avoids overflow or underflow.
softmax <- function (x) {
  y <- exp(x - max(x))
  return(y/sum(y))
}

# Returns the log-density of the normal distribution with zero mean
# and variance s at x.
ldnorm <- function (x, s)
  dnorm(x,sd = sqrt(s),log = TRUE)
