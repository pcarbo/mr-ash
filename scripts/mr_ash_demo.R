# TO DO: Explain here what this script does, and how to use it.

# SIMULATE DATA
# -------------
set.seed(1)
n    <- 500
p    <- 1000
X    <- matrix(rnorm(n*p),n,p)
X    <- scale(X,center = TRUE,scale = FALSE)
sd   <- c(0,0.2,0.5)
w    <- c(0.99,0.05,0.05)
k    <- sample(length(w),p,replace = TRUE,prob = w)
beta <- sd[k] * rnorm(p)
y    <- drop(X %*% beta + rnorm(n))
y    <- y - mean(y)

# FIT MR-ASH MODEL
# ----------------
