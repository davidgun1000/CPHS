library(stochvol)
library(coda)

## Not run:
## Simulate an SV process with leverage
sim <- svsim(10000, mu = -10, phi = 0.95, sigma = 0.2, rho=-0.5)
# Obtain 8000 draws from the sampler (that's too little!)

draws <-
  svsample(sim, draws = 5000, burnin = 100,
           priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.2)



y = read.csv("datareal.csv", header = TRUE)

draws <- svsample(y$y, draws = 10000, burnin = 5000, priormu = c(0, 10000),
                  priorphi = c(100, 1.5), priorsigma = 1, priorrho = c(1, 1))

mu <- draws$para[[1]][,1]
phi <- draws$para[[1]][,2]
tau2 <- (draws$para[[1]][,3])^2
rho <- draws$para[[1]][,5]

h100<-draws$latent[[1]][,100]
IACT_h <- vector(length = 3001)
for (i in 1:3001){
  A=effectiveSize(draws$latent[[1]][,i])
  IACT_h[i] <- length(draws$latent[[1]][,i])/A
  
}

A=effectiveSize(h100)
IACT_h100=length(h100)/A




A=effectiveSize(mu)
IACT_mu=length(mu)/A

A=effectiveSize(phi)
IACT_phi=length(phi)/A

A=effectiveSize(tau2)
IACT_tau2 = length(tau2)/A

A=effectiveSize(rho)
IACT_rho = length(rho)/A
df = data.frame(IACT_mu,IACT_phi,IACT_tau2,IACT_rho)
write.csv(df,"IACT_stochvolrealdata.csv")


summary(draws)
plot(draws)

advanced_draws <-
  svsample(sim, draws = 10000, burnin = 5000,
           priorspec = specify_priors(mu = sv_normal(-10, 1),
                                      sigma2 = sv_gamma(0.5, 2),
                                      rho = sv_beta(4, 4),
                                      nu = sv_constant(5)),
           parallel = "snow", n_chains = 2, n_cpus = 2)

# Example 4
## Predicting USD based on JPY and GBP in the mean
data(exrates)
len <- 3000
ahead <- 30
## Calculate log-returns
logreturns <- apply(exrates[, c("USD", "JPY", "GBP")], 2,
                    function (x) diff(log(x)))
logretUSD <- logreturns[2:(len+1), "USD"]
regressors <- cbind(1, as.matrix(logreturns[1:len, ])) # lagged by 1 day
## Fit SV model to EUR-USD exchange rates
res <- svsample(logretUSD, designmatrix = regressors)