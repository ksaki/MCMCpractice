# Logistic regression with MH and Gibbs

set.seed(123)
library(MCMCpack)

source("logitGibbs.R")
source("logitMH.R")

# simulate data
N <- 1000
X <- runif(N, -2, 2) 
X <- cbind(1, X)
K <- 2
b.true <- runif(K, -2, 2) 

theta <- 1/(1+exp(-X %*% b.true))
Y <- rbinom(N, size=1, prob=theta)

# MLE 
glm.out <- glm(Y ~ -1 + X, family=binomial("logit"))
b.mle <- glm.out$coefficient
b.mle.ci <- confint(glm.out)

# MCMCpack 
iter <- 1000
burnin <- 1000
MCpack.out <- MCMClogit(Y ~ -1 + X,
                        mcmc=iter, burnin=burnin,
                        tune=1)

b.MCpack <- colMeans(MCpack.out)
b.MCpack.ci <- HPDinterval(MCpack.out)

# MH
MH.out <- logitMH(Y=Y, X=X,
                  iter=iter, burnin=burnin)

b.MH <- colMeans(MH.out)
b.MH.ci <- HPDinterval(MH.out)

# Gibbs
gibbs.out <- logitGibbs(Y=Y, X=X,
                        iter=iter, burnin=burnin)

b.gibbs <- colMeans(gibbs.out)
b.gibbs.ci <- HPDinterval(gibbs.out)

# Compare results
cat("Point estimates\n")
cat("true      :", round(b.true, 3),  "\n")
cat("MLE       :", round(b.mle, 3), "\n")
cat("MCMCpack  :", round(b.MCpack, 3), "\n")
cat("MH        :", round(b.MH, 3), "\n")
cat("Gibbs     :", round(b.gibbs, 3), "\n")

cat("Credible/Confidence Interval\n")
print(round(b.mle.ci, 3))
print(round(b.MCpack.ci, 3))
print(round(b.MH.ci, 3))
print(round(b.gibbs.ci, 3))


