# test ordinal IRT 
source("simdata.R")
source("ordinalIRT.R")
set.seed(122)

# Simulation data

# Hyper parameters
N <- 2000 # number of respondents
M <- 10   # number of question itmes
H <- 3   # number of response categories (common to all items)
K <- 1   # number of treatment arms (other than control)

alpha_true <- rep(c(-1, -0.5), 5) 
beta_true <- rep(c(-0.5, -1), 5) 
gamma_true <- matrix(seq(0, H-2), nrow=M, ncol=H-1, byrow=T)
gamma_true <- cbind(-Inf, gamma_true, Inf)

# Generate simulation data
simdata <- get_simdata(N=N,
                       M=M,
                       H=H,
                       K=K,
                       alpha=alpha_true,
                       beta=beta_true,
                       gamma=gamma_true
                       )
print("Response")
print(apply(simdata$Y_mat, 2, table))

# Fit MCMC
out <- main(x=simdata$Y_mat,
            omega.alpha=1,
            omega.beta=1,
            omega.theta=1,  
            omega.gamma=1,
            alpha.init=alpha_true,
            beta.init=beta_true,
            gamma.init=gamma_true,
            theta.init=simdata$theta,
            burnin=0, mcmc=100,
            verbose=10,
            sample.ab=T,
            sample.gamma=T,
            sample.theta=T
            )

# Check results
# gamma
gamma_est <- colMeans(out$gamma[,,-c(1,2,H+1)])
print("gamma")
print(round(gamma_true[,-c(1,2,H+1)], 3))
print(round(gamma_est, 3))

# alpha beta
alpha_est <- colMeans(out$alpha)
beta_est <- colMeans(out$beta)
print("alpha")
print(round(alpha_true, 3))
print(round(alpha_est, 3))
print("beta")
print(round(beta_true, 3))
print(round(beta_est, 3))

# theta
theta_est <- colMeans(out$theta)
plot(simdata$theta, theta_est)

