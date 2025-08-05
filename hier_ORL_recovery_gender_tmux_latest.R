#!/usr/bin/env Rscript
# hier_ORL_recovery_gender_tmux.R
# 1. SETUP: LOAD PACKAGES AND GET SEED

# Get seed from the command line argument passed by the tmux script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  seed <- as.integer(args[1])
} else {
  seed <- 69 # Default seed if not run from the shell script
}
cat(paste("--- Starting Recovery for Seed:", seed, "---\n"))

# Load required packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(R2jags, parallel, extraDistr, truncnorm, tibble, readr)

# Define a function to find the Maximum Posterior Density
MPD <- function(x) {
  density(x)$x[which.max(density(x)$y)]
}

# 2. PAYOFF STRUCTURE: CREATE OR LOAD

# This ensures all parallel tmux panes use the exact same IGT environment.

payoff_file <- "payoff_structure.rds"

if (file.exists(payoff_file)) {
  payoff <- readRDS(payoff_file)
  cat("Loaded existing payoff structure from:", payoff_file, "\n")
} else {
  cat("No payoff structure found. Creating new one with fixed seed (69).\n")
  set.seed(69) # Use a single, fixed seed for generating the task structure
  
  ntrials <- 100
  nstruct <- 10
  
  # Define deck properties
  A_R <- rep(100, nstruct)
  A_L <- c(rep(-250, nstruct * 0.5), rep(0, nstruct * 0.5))
  B_R <- rep(100, nstruct)
  B_L <- c(rep(-1250, nstruct * 0.1), rep(0, nstruct * 0.9))
  C_R <- rep(50, nstruct)
  C_L <- c(rep(-50, nstruct * 0.5), rep(0, nstruct * 0.5))
  D_R <- rep(50, nstruct)
  D_L <- c(rep(-250, nstruct * 0.1), rep(0, nstruct * 0.9))
  
  # Create pseudorandomized full payoff structure
  A <- B <- C <- D <- array(NA, ntrials)
  for (i in 1:(ntrials / nstruct)) {
    idx <- (1 + (i - 1) * nstruct):(i * nstruct)
    A[idx] <- A_R + sample(A_L)
    B[idx] <- B_R + sample(B_L)
    C[idx] <- C_R + sample(C_L)
    D[idx] <- D_R + sample(D_L)
  }
  
  payoff <- cbind(A, B, C, D) / 100
  saveRDS(payoff, file = payoff_file)
  cat("Saved new payoff structure to:", payoff_file, "\n")
}

# Now, set the seed for this specific parallel run
set.seed(seed)

# 3. RECOVERY SIMULATION

# --- Simulation Settings ---
n_iterations <- 24 # Total iterations for this seed
n_subs <- 139      # 98 males, 41 females
ntrials_all <- rep(100, n_subs)

# --- Data Collection ---
# Initialize  empty tibble to store results iter by iter
output_df <- tibble()

cat("---- Starting recovery loop ----\n")
start_time <- Sys.time()
print(start_time)

for (i in 1:n_iterations) {
  
  # --- Generate True Parameters for this Iteration ---
  # Group-level means (mu)
  mu_a_rew   <- runif(1, 0, 1)
  mu_a_pun   <- runif(1, 0, 1)
  mu_K       <- runif(1, 0, 2)
  mu_omega_f <- runif(1, -2, 2)
  mu_omega_p <- runif(1, -2, 2)
  
  # Group-level standard deviations (sigma)
  sigma_a_rew   <- runif(1, 0.1, 0.5)
  sigma_a_pun   <- runif(1, 0.1, 0.5)
  sigma_K       <- runif(1, 0.1, 0.5)
  sigma_omega_f <- runif(1, 0.1, 0.5)
  sigma_omega_p <- runif(1, 0.1, 0.5)
  
  # Gender difference effects (alpha)
  alpha_a_rew   <- runif(1, -1, 1) # Centered on 0
  alpha_a_pun   <- runif(1, -1, 1) # Centered on 0
  alpha_omega_f <- runif(1, -2, 2) # Centered on 0
  
  # --- Simulate IGT Data ---
  source('hier_ORL_gender_sim_latest.R')
  ORL_sims <- hier_ORL_gender_sim(payoff, n_subs, ntrials_all,
                                  mu_a_rew, alpha_a_rew,
                                  mu_a_pun, alpha_a_pun,
                                  mu_K,
                                  mu_omega_f, alpha_omega_f,
                                  mu_omega_p,
                                  sigma_a_rew, sigma_a_pun, sigma_K,
                                  sigma_omega_f, sigma_omega_p)
  
  x      <- ORL_sims$x
  X      <- ORL_sims$X
  gender <- ORL_sims$gender
  
  # --- Fit Model with JAGS ---
  jags_data <- list(x = x, 
                    X = X, 
                    ntrials = ntrials_all, 
                    nsubs = n_subs, 
                    gender = gender)
  
  jags_params <- c("mu_a_rew", "alpha_a_rew", "mu_a_pun", "alpha_a_pun",
                   "mu_K", "mu_omega_f", "alpha_omega_f", "mu_omega_p",
                   "lambda_a_rew", "lambda_a_pun", "lambda_K",
                   "lambda_omega_f", "lambda_omega_p")
  
  samples <- jags.parallel(data = jags_data, inits = NULL, parameters.to.save = jags_params,
                           model.file = "hier_ORL_gender_latest.txt",
                           n.chains = 4, n.iter = 5000, n.burnin = 1000,
                           n.thin = 1, n.cluster = 4)
  
  # Extract posterior simulations
  Y <- samples$BUGSoutput$sims.list
  
  # --- Store Results for this Iteration ---
  temp_df <- tibble(
    # General Info
    seed = seed,
    iteration = i,
    
    # True Mu Parameters
    true_mu_a_rew = mu_a_rew,
    true_mu_a_pun = mu_a_pun,
    true_mu_K = mu_K,
    true_mu_omega_f = mu_omega_f,
    true_mu_omega_p = mu_omega_p,
    
    # Inferred Mu Parameters (MPD)
    infer_mu_a_rew = MPD(Y$mu_a_rew),
    infer_mu_a_pun = MPD(Y$mu_a_pun),
    infer_mu_K = MPD(Y$mu_K),
    infer_mu_omega_f = MPD(Y$mu_omega_f),
    infer_mu_omega_p = MPD(Y$mu_omega_p),
    
    # True Alpha (Gender Effect) Parameters
    true_alpha_a_rew = alpha_a_rew,
    true_alpha_a_pun = alpha_a_pun,
    true_alpha_omega_f = alpha_omega_f,
    
    # Inferred Alpha Parameters (MPD)
    infer_alpha_a_rew = MPD(Y$alpha_a_rew),
    infer_alpha_a_pun = MPD(Y$alpha_a_pun),
    infer_alpha_omega_f = MPD(Y$alpha_omega_f),
    
    # True Sigma Parameters (for comparison)
    true_sigma_a_rew = sigma_a_rew,
    true_sigma_a_pun = sigma_a_pun,
    true_sigma_K = sigma_K,
    true_sigma_omega_f = sigma_omega_f,
    true_sigma_omega_p = sigma_omega_p,
    
    # Inferred Sigma Parameters (derived from lambda)
    infer_sigma_a_rew = 1/sqrt(MPD(Y$lambda_a_rew)),
    infer_sigma_a_pun = 1/sqrt(MPD(Y$lambda_a_pun)),
    infer_sigma_K = 1/sqrt(MPD(Y$lambda_K)),
    infer_sigma_omega_f = 1/sqrt(MPD(Y$lambda_omega_f)),
    infer_sigma_omega_p = 1/sqrt(MPD(Y$lambda_omega_p))
  )
  
  # Append the results of the current iteration to the main data frame
  output_df <- rbind(output_df, temp_df)
  
  # --- Progress Update ---
  cat(sprintf("Seed %d: Progress %.1f%% complete. Time elapsed: %s\n", 
              seed, (i / n_iterations * 100), format(Sys.time() - start_time)))
}

# 4. SAVE FINAL RESULTS

# Create a unique filename for this seed's results and save the data frame.

# Create a 'results' directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

output_filename <- paste0("results/recovery_results_seed_", seed, ".rds")
write_rds(output_df, output_filename)

cat(paste("\n--- Finished recovery for seed", seed, "---"))
cat(paste("\nTotal time:", format(Sys.time() - start_time)))
cat(paste("\nResults saved to:", output_filename, "\n"))