# functions for logistic regression with MH
library(mvtnorm)

# prior
log.prior <- function(b, b0, B0){
  K <- length(b)
  out <- dmvnorm(b, mean=b0, sigma=solve(B0), log=T) 
  return(out)
}

# likelihood
log.like <- function(b, Y, X){
  theta <- 1/(1+exp(-X %*% b))
  out <- dbinom(Y, size=1, prob=theta, log=TRUE)
  return(sum(out))
}

# posterior
log.post.kernel <- function(b, Y, X, b0, B0){
  return(log.like(b, Y, X) + log.prior(b, b0, B0))
}

# Logistic Regression MH with diagonal covariance matrix 
# for proposal density
logitMH <- function(Y, # outcome
                    X, # covariate (incluidng intercept)
                    b0=rep(0, ncol(X)), # prior mean
                    B0=diag(0.001, ncol(X)), # prior precision 
                    b.init=rep(0, ncol(X)), # initial values
                    iter=1000,
                    burnin=1000,
                    sigma.cand=diag(1, ncol(X))){

  total.iter <- iter + burnin
  K <- ncol(X)

  # setup storage matrix for theta draws
  b.store <- matrix(NA, total.iter, K)

  # initial value of theta 
  b <- b.init

  # acceptance counter
  accepts <- 0

  # the Metropolis sampling begins here
  for (m in seq(1, total.iter)){

    # generate candidate
    # Use hessian for covariance matrix
    b.can <- b + as.numeric(rmvnorm(1, mean=rep(0,K),
                                    sigma=sigma.cand))

    # calculate acceptance ratio
    ratio <- exp(log.post.kernel(b.can, Y, X, b0, B0) -
                 log.post.kernel(b, Y, X, b0, B0))
    alpha <- min(ratio, 1)

    # accept or reject candidate based on alpha
    if (runif(1) < alpha){
      b <- b.can
      accepts <- accepts + 1
    }
    
    # store sample
    b.store[m,] <- b

    # counter
    if (m %% 1000 == 0){
      cat("iteration: ", m, "\n")
      cat("Acceptance Ratio: ", round(accepts/m, 3), "\n")
    }
  }
  return(mcmc(b.store))
}



# Logistic Regression MH with Hessian to control covariance of 
# proposal density
logitMH.H <- function(Y, # outcome
                      X, # covariate (inclidng intercept)
                      b0=rep(0, ncol(X)), # prior mean
                      B0=diag(0.001, ncol(X)), # prior precision 
                      b.init=rep(0, ncol(X)), # initial values
                      iter=1000,
                      burnin=1000,
                      tune=1.1){

  total.iter <- iter + burnin
  K <- ncol(X)

  # setup storage matrix for theta draws
  b.store <- matrix(NA, total.iter, K)

  # initial value of theta 
  b <- b.init

  # acceptance counter
  accepts <- 0

  # Use hessian (variance-covariance matrix from MLE)
  # to give direction in proposal distribution
  # Posterior variance with some tuning 
  glm.out <- glm(Y ~ -1 + X, family=poisson("log"))
  V <- vcov(glm.out)
  propV <- tune * solve(B0 + solve(V)) * tune

  # the Metropolis sampling begins here
  for (m in seq(1, total.iter)){

    # generate candidate
    # Use hessian for covariance matrix
    b.can <- b + as.numeric(rmvnorm(1, mean=rep(0,K),
                                    sigma=propV))

    # calculate acceptance ratio
    ratio <- exp(log.post.kernel(b.can, Y, X, b0, B0) -
                 log.post.kernel(b, Y, X, b0, B0))
    alpha <- min(ratio, 1)

    # accept or reject candidate based on alpha
    if (runif(1) < alpha){
      b <- b.can
      accepts <- accepts + 1
    }
    
    # store sample
    b.store[m,] <- b

    # counter
    if (m %% 1000 == 0){
      cat("iteration: ", m, "\n")
      cat("Acceptance Ratio: ", round(accepts/m, 3), "\n")
    }
  }
  return(mcmc(b.store))
}


