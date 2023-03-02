# functions for logistic regression Gibbs sampler
library(BayesLogit)
#library(pgdraw) # alternative  Polya-Gamma sampler

# Logistic Regression Gibbs sampler with Polya-Gamma auxillary variable
logitGibbs <- function(Y, # outcome
                    X, # covariate (incluidng intercept)
                    b0=rep(0, ncol(X)), # prior mean
                    B0=diag(100000, ncol(X)), # prior variance
                    b.init=NULL, # initial values
                    omega.init=diag(rep(1, nrow(X))),
                    iter=1000,
                    burnin=1000,
                    verbose=100
                    ){

  total.iter <- iter + burnin
  K <- ncol(X)
  N <- nrow(X)
  B0.inv <- solve(B0)
  k <- Y-1/2

  # setup storage matrix for theta draws
  b.store <- matrix(NA, iter, K)

  # initial value
  glm.fit <- glm(Y ~ -1 + X, family=binomial("logit"))
  b <- glm.fit$coefficients
  omega <- omega.init

  # Gibbs 
  for (m in seq(1, total.iter)){

    # sample beta given omega
    omega <- diag(rpg(N, rep(1, N), X %*% b))

    # sample omega given beta
    V.beta <- solve(t(X) %*% omega %*% X + B0.inv)
    m.beta <- V.beta %*% (crossprod(X, k) + B0.inv %*% b0)
    b <- as.numeric(rmvnorm(1, m.beta, V.beta))
    
    
    # store sample
    if (m > burnin){
      b.store[m-burnin,] <- b
    }

    # counter
    if (m %% verbose == 0){
      cat("iteration: ", m, " / ", total.iter, "\n")
    }
  }
  return(mcmc(b.store))
}



