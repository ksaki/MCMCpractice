# MCMC for Ordinal IRT

library(mvtnorm)
library(truncnorm)
main <- function(x, 
                 omega.alpha=1,
                 omega.beta=1,
                 omega.theta=1,
                 omega.gamma=1,
                 alpha.init=NULL,
                 beta.init=NULL,
                 gamma.init=NULL,
                 theta.init=NULL,
                 burnin=100,
                 mcmc=1000,
                 verbose=100,
                 tune=1,
                 sample.ab=T,
                 sample.gamma=T,
                 sample.theta=T
                 ){

    # Note:
    # gamma is M by H+1 matrix. The first col is -Inf, 
    # the last col is Inf. The second col is 1.

    # utility objects
    total.iter <- burnin + mcmc
    N <- nrow(x) # observation
    M <- ncol(x) # item
    H <- max(x)  # response category
    accept <- rep(0, M)  # acceptance counter

    # storage
    out.alpha <- matrix(NA, total.iter, M)
    out.beta <- matrix(NA, total.iter, M)
    out.gamma <- array(NA, dim=c(total.iter, M, H+1))
    out.theta <- matrix(NA, total.iter, N)

    # initailiation
    alpha <- alpha.init
    beta <- beta.init
    gamma <- gamma.init
    theta <- theta.init
    xstar <- matrix(0, nrow=N, ncol=M)
    xstar.mean <- matrix(0, nrow=N, ncol=M)

    # mcmc loop
    for (iter in seq(1, total.iter)){
      

      ################################################################
      ## sample theta
      if (sample.theta){
        for (i in seq(1, N)){
          xstar.theta <- xstar[i,] + alpha 
          sigma.theta <- solve(crossprod(beta) + 1/omega.theta)
          m.theta <- sigma.theta * crossprod(beta,  xstar.theta) 
          theta[i] <- rnorm(1, m.theta, sqrt(sigma.theta))
        }
      }

      ################################################################
      # sample alpha and beta
      if (sample.ab){
        ab <- matrix(NA, nrow=M, ncol=2)
        theta.ab <- cbind(-1, theta)
        theta.ab.cp <- crossprod(theta.ab)
        omega.ab.inv <- solve(diag(c(omega.alpha, omega.beta)))
        sigma.ab <- solve(theta.ab.cp + omega.ab.inv)
        for (j in seq(1, M)){
          m.ab <- sigma.ab %*% crossprod(theta.ab, xstar[,j])
          ab[j,] <- rmvnorm(1, mean=m.ab, sigma=sigma.ab)

          # all beta is constrained to be positive
          while (any(ab[j,2] > 0)){
            ab[j,] <- rmvnorm(1, mean=m.ab, sigma=sigma.ab)
          }

        }
        alpha <- ab[,1]
        beta <- ab[,2]
      }

      ###############################################################
      # sample xstar 
      for (i in seq(1, N)){
        for (j in seq(1, M)){
          xstar.mean[i,j] <- - alpha[j] + beta[j] * theta[i]
          xstar[i,j] <- rtruncnorm(n=1, mean=xstar.mean[i,j], sd=1,
                                   a=gamma[j,x[i,j]],
                                   b=gamma[j,x[i,j]+1]) 
        }
      }

      ################################################################
      # sample gamma (MH)
      if (H == 2){
        # no gamma update required if binary responses

      } else if (H > 2){
        for (j in seq(1, M)){
          gamma.prop <- rep(0, H+1)
          gamma.prop[1] <- -Inf
          gamma.prop[2] <- 0
          gamma.prop[H+1] <- Inf
          for (h in seq(3, H)){
            gamma.prop[h] <- rtruncnorm(1,
                                        a=gamma.prop[h-1],
                                        b=gamma[j,h+1],
                                        mean=gamma[j,h],
                                        sd=tune)
          }
        
          # acceptance ratio
          rho <- 0
          prop.ratio.log <- 0
          like.ratio.log <- 0
          for (h in seq(2, H-1)){
            prop.ratio.log <- prop.ratio.log + 
              log(pnorm((gamma[j,h+2] - gamma[j,h+1])/omega.gamma) -
                  pnorm((gamma.prop[h] - gamma[j,h+1])/omega.gamma)) -
              log(pnorm((gamma.prop[h+2] - gamma.prop[h+1])/omega.gamma)-
                 pnorm((gamma[j,h] - gamma.prop[h+1])/omega.gamma))
          }

          for (i in seq(1, N)){
            like.ratio.log <- like.ratio.log +
              log(pnorm(gamma.prop[x[i,j]+1] - xstar.mean[i,j]) -
                  pnorm(gamma.prop[x[i,j]] - xstar.mean[i,j])) -
              log(pnorm(gamma[j,x[i,j]+1] - xstar.mean[i,j]) -
                  pnorm(gamma[j,x[i,j]] - xstar.mean[i,j]))
          }
          rho <- exp(prop.ratio.log + like.ratio.log)
          if (runif(n=1, min=0, max=1) < rho){
            gamma[j,] <- gamma.prop
            accept[j] <- accept[j] + 1
          }
        }
      }

      # store results
      out.alpha[iter,] <- alpha
      out.beta[iter,] <- beta
      out.gamma[iter,,] <- gamma
      out.theta[iter,] <- theta

      if (iter %% verbose == 0){
        cat("iteration = ", iter, " / ", total.iter, "\n")
        cat("acceptance ratio:\n",
            round(accept/iter, 2), "\n")
      }
    }
    return(list(alpha=out.alpha,
                beta=out.beta,
                gamma=out.gamma,
                theta=out.theta
                ))

}
