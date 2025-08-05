# CONVERGENCE DIAGNOSTICS
install.packages('coda')
library(coda)

# Load results file
setwd('/work/LinusJoakimBackström#7558/DM exam')
load("results/ORL_gender_latest_results.RData") # This loads the 'fit' object

print(fit)
# highest Rhat value is 1.004

# Build mcmc.list manually from sims.array
samples_array <- fit$BUGSoutput$sims.array  # [iterations, chains, parameters]
param_names <- dimnames(samples_array)[[3]]

# Convert each chain to mcmc
mcmc_chains <- lapply(1:dim(samples_array)[2], function(chain) {
  mcmc(samples_array[, chain, ])
})

mcmc_samples <- mcmc.list(mcmc_chains)


plot(mcmc_samples[, c("alpha_a_pun", "alpha_a_rew", "alpha_omega_f")])