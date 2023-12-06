library("rstan")
#T=3000
#mu = -1
#phi = 0.98
#tau2 = 0.1
#rho = -0.5

#ht <- matrix(NA,T,1)
#y <- matrix(NA,T,1)

#ht[1] <- sqrt(tau2/(1-phi^2))*rnorm(1)
#y[1] <- sqrt(exp(ht[1]))*rnorm(1)

#for(t in 2:T){
#   ht[t] <- mu + phi*(ht[t-1]-mu) + rho*sqrt(tau2)*exp(-ht[t-1]/2)*y[t-1] + sqrt(tau2)*sqrt(1-rho^2)*rnorm(1)  
#   y[t] <- sqrt(exp(ht[t]))*rnorm(1)
#}

T=3001
y = read.csv("datareal.csv", header = TRUE)
compiled_model <- stan_model("stochasticvolatilityLeverage.stan")
data_list <- list(T = T, y = y)

fitted_model <- sampling(compiled_model, data = data_list, cores = 1, chains = 1, iter = 15000, control=list(max_treedepth=10,stepsize = 1))
hmc.fit <- extract(fitted_model, pars = c("mu","phi","sigma","h"),
                   permuted = F)

plot(hmc.fit[,,1],type = "l")
plot(hmc.fit[,,2],type = "l")
plot(hmc.fit[,,3],type = "l")



#N <- 100L #number of individuals
#n <- 10L # number of responses per individual
#G <- N*n
#PX <- 1 # Number of exogenous variables
#PN <- 1 # Number of endogenous variables
#PZ <- 2 # Number of instruments
#compiled_model <- stan_model("linear_fixedeffect_endo.stan")
# Exogenous variables (make them correlated with a random correlation matrix)
#X_exog <- MASS::mvrnorm(N, rep(0, PX), diag(PX))

# Some fake instruments
#Z <- MASS::mvrnorm(N, rep(0, PZ), diag(PZ))

#param_A <- 0.1
#param_gamma <- 0.1
#param_delta_endo <- c(0.1,0.1)
#param_delta_exo <- 0.1
#param_beta <- 0.1
#param_a <- 0.1
#n_fixed_effects <- length(c())
#tau <- 0.9
#omega <- matrix(c(2,0.2,0.2,2),2,2)
#X_exog <- list()
#X_endo <- list()
#Z <- list()

#y_outcome <- list()

#for (i in 1:N) {
#    X_exog[[i]] <- matrix(rnorm(n * PX), nrow = n, ncol = PX)
#    alpha_i <- rnorm(1, 0, 1)
#    Z[[i]] <- matrix(rnorm(n * PZ), nrow = n, ncol = PZ)
#    
#    err <- MASS::mvrnorm(n, rep(0, (1+PN)), omega)
#    X_endo[[i]] <- param_A + X_exog[[i]] %*% param_gamma +  Z[[i]] %*% param_delta_endo + err[,2]
     
#    y_outcome[[i]] <- alpha_i + param_a + X_endo[[i]] %*% param_beta + X_exog[[i]] %*% param_delta_exo + err[,1]
    
    
    #p_i <- exp(X[[i]] %*% beta + alpha_i) / (1 + exp(X[[i]] %*% beta + alpha_i))
    #y[[i]] <- rbinom(n, 1, p_i)
#}


  


## Data manipulation ##
#y_long <- unlist(y_outcome) #as.vector(t(y))
#X_exog_long <- do.call("rbind",X_exog)
#X_endo_long <- do.call("rbind",X_endo)
#Z_long <- do.call("rbind",Z)
#group <- rep(1:N, each = n)

#data_list <- list(G = G, N = N, PX = PX, PN = PN, PZ = PZ, 
#                 X_exog = X_exog_long, X_endog = X_endo_long, 
#                 Z = Z_long, Y_outcome = y_long, g = group)

#fitted_model <- sampling(compiled_model, data = data_list, cores = 1, chains = 1, iter = 1000, control=list(max_treedepth=10,stepsize = 1))

#hmc.fit <- extract(fitted_model, pars = c("gamma1","gamma2","alpha","L_Omega","scale"),
#                   permuted = F)

#hmc.fit.L_Omega <- extract(fitted_model, pars = c("L_Omega"), permuted = F)
#hmc.fit.scale <- extract(fitted_model, pars = c("scale"), permuted = F)

#plot(hmc.fit[,,1],type = "l")
#plot(hmc.fit[,,2],type = "l")
#plot(hmc.fit[,,3],type = "l")
#plot(hmc.fit[,,4],type = "l")
#plot(hmc.fit[,,5],type = "l")
#plot(hmc.fit[,,6],type = "l")
#plot(hmc.fit[,,7],type = "l")
#plot(hmc.fit[,,8],type = "l")
#plot(hmc.fit[,,9],type = "l")
#plot(hmc.fit[,,10],type = "l")
#plot(hmc.fit[,,11],type = "l")

#hmc.samples.L_Omega <- matrix(NA, 1000, 4)
#for (p in 1:4) {
#  hmc.samples.L_Omega[, p] <- hmc.fit.L_Omega[, , p]
#}
#hmc.post_mean.L_Omega <- as.vector(apply(hmc.samples.L_Omega, 2, mean))

#hmc.samples.scale <- matrix(NA, 1000, 2)
#for (p in 1:2) {
#  hmc.samples.scale[, p] <- hmc.fit.scale[, , p]
#}
#hmc.post_mean.scale <- as.vector(apply(hmc.samples.scale, 2, mean))


#hfit <- run_stan_logmm(iters = hmc.iters, data = y_long, 
#                       grouping = rep(1:N, each = n), n_groups = N,
#                       fixed_covariates = X_long)

#if (save_hmc_results) {
#  saveRDS(hfit, file = paste0("logistic_mm_hmc_N", N, "_n", n, "_", date, ".rds"))
#}
 