# Generate simulation data for ordinal IRT

get_simdata <- function(N=2000,
                        M=2,
                        K=1,
                        H=5,
                        alpha=NULL,
                        beta=NULL,
                        gamma=NULL,
                        theta=NULL,
                        omega_alpha=1,
                        omega_beta=1,
                        omega_tau=1
                        ){

  # N (int) number of respondents
  # M (int) number of policies
  # K (int) number of treatment types (excluding control) 
  # H (int) response category
  # alpha (vector) M difficulty parameters
  # beta (vector) M discrimination parameters

  # Draw theta
  if (is.null(theta)){
    theta_vec <- rnorm(N, 0, 1)
    theta_vec <- as.matrix(theta_vec, nrow=N, ncol=1) 
  } else {
    theta_vec <- theta
  }

  # Draw difficulty
  if (is.null(alpha)){
    alpha_vec <- rnorm(n=M, mean=0, sd=sqrt(omega_alpha))
  } else {
    alpha_vec <- alpha
  }
  
  # Draw discrimination
  if (is.null(beta)){
    beta_vec <- rnorm(n=M, mean=0, sd=sqrt(omega_beta))
  } else {
    beta_vec <- beta
  }
  
  # Draw cutoff
  if (is.null(gamma)){
    gamma_mat <- matrix(0, M, H-1)
    for (j in seq(1, M)){
      gamma_mat[j,] <- sort(runif(H-1, -3, 3))
    }
  } else {
    gamma_mat <- gamma 
  }

  # Draw resopnses
  Y_mat <- matrix(NA, nrow=N, ncol=M) # response matrix
  Y_star_mat <- matrix(NA, nrow=N, ncol=M) 
  for (i in seq(1, N)){
    for (j in seq(1, M)){
      # Get the linear part based on treatment status
      ln <- - alpha_vec[j] + beta_vec[j] * theta_vec[i] 

      # auxillary variable
      Y_star_mat[i,j] <- rnorm(1, mean=ln, sd=1)

      # for ordered case
      for (h in seq(1, H)){
        if (gamma_mat[j,h] < Y_star_mat[i,j] & Y_star_mat[i,j] < gamma_mat[j,h+1]){
          Y_mat[i,j] <- h
          break
        }
      }
    }
  }

  return(list(Y_mat=Y_mat, 
              Y_star_mat=Y_star_mat,
              alpha_vec=alpha_vec,
              beta_vec=beta_vec,
              gamma_mat=gamma_mat,
              theta_vec=theta_vec,
              N=N,
              M=M,
              K=K,
              H=H,
              omega_alpha=omega_alpha,
              omega_beta=omega_beta
              ))
}



