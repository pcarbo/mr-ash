# TO DO: Explain here what this script is for.

# Script settings.
n  <- 75
s  <- 3
b  <- c(-1,1,0,0)
s0 <- c(0.001,0.01,0.1,1)
w0 <- c(0.0001,0.9,0.05,0.05)
w0 <- w0/sum(w0)

# Simulate data.
set.seed(1)
p <- length(b)
X <- matrix(rnorm(n*p),n,p)
X <- scale(X,center = TRUE,scale = FALSE)
y <- drop(X %*% b) + sqrt(s)*rnorm(n)
y <- y - mean(y)

# Fit model using varbvsmix.
s00    <- s0
s00[1] <- 0
fit <- varbvsmix(X,Z = NULL,y,(s00/s)^2,s^2,w0,update.sigma = FALSE,
                 update.sa = FALSE,update.w = FALSE,verbose = TRUE,
                 tol = 1e-8)
fit$logZ <- fit$logZ + log(n)/2
best <- rowSums(fit$alpha * fit$mu)

# Fit model using optim.
norm2 <- function (x)
  sqrt(sum(x^2))
betavarmix <- function (p, mu, s)
  rowSums(p*(s + mu^2)) - rowSums(p*mu)^2
normalizelogweights <- function (logw) {
  c <- max(logw)
  w <- exp(logw - c)
  return(w/sum(w))
}
compute_elbo <- function (bhat, X, y, sigma, s0, w) {
  n     <- length(y)
  p     <- length(bhat)
  K     <- length(s0)
  d     <- diag(crossprod(X))
  alpha <- matrix(0,p,K)
  mu    <- matrix(0,p,K)
  s     <- matrix(0,p,K)
  for (i in 1:p) {
    shat   <- sigma/sqrt(d[i])
    s1     <- 1/sqrt(1/shat^2 + 1/s0^2)
    mu1    <- (s1/shat)^2*bhat[i]
    logBF  <- log(shat/sqrt(shat^2 + s0^2)) + (mu1/s1)^2/2
    logp1  <- log(w) + logBF
    mu[i,] <- mu1
    s[i,]  <- s1^2
    alpha[i,] <- normalizelogweights(logp1)
  }
  b <- rowSums(alpha * mu)
  print(b)
  elbo <- (-n/2*log(2*pi*sigma^2)
           - (norm2(y - X %*% b)^2 
              + sum(d*betavarmix(alpha,mu,s)))/(2*sigma^2))
  for (i in 1:K)
    elbo <- (elbo + sum(alpha[,i]*log(w[i] + 1e-15)) 
                  - sum(alpha[,i]*log(alpha[,i] + 1e-15)))
  for (i in 1:K)
    elbo <- (elbo+(sum(alpha[,i]) + sum(alpha[,i]*log(s[,i]/(s0[i]^2))))/2
             - sum(alpha[,i]*(s[,i] + mu[,i]^2))/(s0[i]^2)/2)
  print(elbo,digits = 8)
  return(elbo)
}
compute_elbo2 <- function (bhat, X, y, sigma, s0, w) {
  n     <- length(y)
  p     <- length(bhat)
  K     <- length(s0)
  d     <- diag(crossprod(X))
  alpha <- matrix(0,p,K)
  mu    <- matrix(0,p,K)
  s     <- matrix(0,p,K)
  shat  <- rep(0,p)
  llik  <- rep(0,p)
  for (i in 1:p) {
    shat[i] <- sigma/sqrt(d[i])
    s1     <- 1/sqrt(1/shat[i]^2 + 1/s0^2)
    mu1    <- (s1/shat[i])^2*bhat[i]
    logBF  <- log(shat[i]/sqrt(shat[i]^2 + s0^2)) + (mu1/s1)^2/2
    logp1  <- log(w) + logBF
    u      <- max(logp1)
    llik[i] <- sum(dnorm(y,0,sigma,log = TRUE)) +
               log(sum(exp(logp1 - u))) + u
    mu[i,] <- mu1
    s[i,]  <- s1^2
    alpha[i,] <- normalizelogweights(logp1)
  }
  b <- rowSums(alpha * mu)
  print(b)
  elbo <- n*(p-1)/2*log(2*pi*sigma^2) - norm2(y - X %*% b)^2/(2*sigma^2) +
          sum(llik) + sum(((bhat - b)^2 + sum(y^2)/d - bhat^2)/(2*shat^2))
  return(elbo)
}
theta0 <- rnorm(p)
out <- optim(theta0,function (par) -compute_elbo2(par,X,y,s,s0,w0),
             method = "Nelder-Mead",control = list(reltol = 1e-10))
print(best)
print(max(fit$logZ),digits = 10)
print(-out$value,digits = 10)
