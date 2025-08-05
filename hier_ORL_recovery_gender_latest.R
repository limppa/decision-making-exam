#!/usr/bin/env Rscript
# hier_ORL_recovery_gender_tmux.R

# -----------------------------
# 1. Get the seed from the command line 
#    (as done in Tmux_ORL_rec.R)
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  seed <- as.numeric(args[1])
} else {
  seed <- 123  # default if none provided
}
cat("Using seed:", seed, "\n")

# -----------------------------
# 2. Load packages and set working directory
# -----------------------------
if (!require(pacman)) install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm, ggplot2)

setwd('/work/LinusJoakimBackström#7558/DM exam')

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

# -----------------------------
# 4. Create (or load) the payoff structure 
#    so that all sub‐runs use the same one.
#    We use a fixed seed (e.g., 123) when generating it.
# -----------------------------
payoff_file <- "payoff_structure.rds"
if (file.exists(payoff_file)) {
  payoff <- readRDS(payoff_file)
  cat("Loaded existing payoff structure from", payoff_file, "\n")
} else {
  cat("No existing payoff structure found. Creating a new one...\n")
  set.seed(69)  # fixed seed for payoff structure generation
  
  #------ Create task environment -------
  ntrials <- 100
  nstruct <- 10
  freq <- 0.5    # Loss frequency for frequent loss decks (A/C)
  infreq <- 0.1  # Loss frequency for infrequent loss decks (B/D)
  
  # Deck A (bad frequent loss): variable losses around mean -250
  bad_freq_l <- c(-150, -200, -250, -300, -350)  # Mean = -250
  bad_r <- 100  # winnings for "bad" decks (A&B)
  
  # Deck B (bad infrequent loss)
  bad_infreq_l <- -1250
  
  # Deck C (good frequent loss): variable losses around mean -50
  good_freq_l <- c(-25, -50, -75)  # Mean = -50
  good_r <- 50   # winnings for "good" decks (C&D)
  
  # Deck D (good infrequent loss)
  good_infreq_l <- -250
  
  # Create loss sequences for each deck type
  A_L <- c(sample(bad_freq_l, nstruct * freq, replace = TRUE),
           rep(0, nstruct * (1 - freq)))
  B_L <- c(rep(bad_infreq_l, nstruct * infreq),
           rep(0, nstruct * (1 - infreq)))
  C_L <- c(sample(good_freq_l, nstruct * freq, replace = TRUE),
           rep(0, nstruct * (1 - freq)))
  D_L <- c(rep(good_infreq_l, nstruct * infreq),
           rep(0, nstruct * (1 - infreq)))
  
  # Build payoff structure
  A <- array(NA, ntrials)
  B <- array(NA, ntrials)
  C <- array(NA, ntrials)
  D <- array(NA, ntrials)
  
  for (i in 1:(ntrials / nstruct)) {
    idx <- (1 + (i - 1) * nstruct):(i * nstruct)
    A[idx] <- rep(bad_r, nstruct) + sample(A_L)
    B[idx] <- rep(bad_r, nstruct) + sample(B_L)
    C[idx] <- rep(good_r, nstruct) + sample(C_L)
    D[idx] <- rep(good_r, nstruct) + sample(D_L)
  }
  
  # Combine decks as columns and rescale (dividing by 100)
  payoff <- cbind(A, B, C, D) / 100
  
  # Save the payoff structure so that all sub-runs use the same one.
  saveRDS(payoff, payoff_file)
  cat("Created and saved new payoff structure to", payoff_file, "\n")
}

# Now set the seed for the parameter recovery (based on command-line argument)
set.seed(seed)

# -----------------------------
# 5. Set up parameter recovery settings (for the gender version)
# -----------------------------
niterations <- 100
nsubs <- 139         # mimicking the data structure from Kildahl et al.
ntrials_all <- rep(100, nsubs)  # each subject has 100 trials

# Preallocate arrays for group-level parameters (mu) ...
true_mu_a_rew   <- array(NA, niterations)
true_mu_a_pun   <- array(NA, niterations)
true_mu_K       <- array(NA, niterations)
true_mu_omega_f <- array(NA, niterations)
true_mu_omega_p <- array(NA, niterations)

infer_mu_a_rew   <- array(NA, niterations)
infer_mu_a_pun   <- array(NA, niterations)
infer_mu_K       <- array(NA, niterations)
infer_mu_omega_f <- array(NA, niterations)
infer_mu_omega_p <- array(NA, niterations)

# ... for precision (lambda, computed from sigma)
true_lambda_a_rew   <- array(NA, niterations)
true_lambda_a_pun   <- array(NA, niterations)
true_lambda_K       <- array(NA, niterations)
true_lambda_omega_f <- array(NA, niterations)
true_lambda_omega_p <- array(NA, niterations)

infer_lambda_a_rew   <- array(NA, niterations)
infer_lambda_a_pun   <- array(NA, niterations)
infer_lambda_K       <- array(NA, niterations)
infer_lambda_omega_f <- array(NA, niterations)
infer_lambda_omega_p <- array(NA, niterations)

# ... and for gender effects (alpha)
true_alpha_a_rew   <- array(NA, niterations)
true_alpha_a_pun   <- array(NA, niterations)
true_alpha_omega_f <- array(NA, niterations)

infer_alpha_a_rew   <- array(NA, niterations)
infer_alpha_a_pun   <- array(NA, niterations)
infer_alpha_omega_f <- array(NA, niterations)

# -----------------------------
# 6. Run the recovery iterations
# -----------------------------
cat("---- First iteration starts ----\n")
start_time <- Sys.time()
cat("Start time:", start_time, "\n\n")

for (i in 1:niterations) {
  cat("Iteration:", i, "\n")
  
  ntrials <- ntrials_all
  
  # Draw group-level (mu) and individual variability (sigma) parameters
  mu_a_rew   <- runif(1, 0, 1)
  mu_a_pun   <- runif(1, 0, 1)
  mu_K       <- runif(1, 0, 2)
  mu_omega_f <- runif(1, -2, 2)
  mu_omega_p <- runif(1, -2, 2)
  
  sigma_a_rew   <- runif(1, 0.1, 0.5)
  sigma_a_pun   <- runif(1, 0.1, 0.5)
  sigma_K       <- runif(1, 0.1, 0.5)
  sigma_omega_f <- runif(1, 0.1, 0.5)
  sigma_omega_p <- runif(1, 0.1, 0.5)
  
  # Draw gender effect parameters
  alpha_a_rew   <- runif(1, 0, 1)
  alpha_a_pun   <- runif(1, 0, 1)
  alpha_omega_f <- runif(1, -2, 2)
  
  # Load and run the simulation function (which uses your gender‐specific ORL simulation)
  source('hier_ORL_gender_sim_latest.R')
  ORL_sims <- hier_ORL_gender_sim(payoff, nsubs, ntrials,
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
  
  # Set up JAGS: data and parameters to monitor
  data_jags <- list("x", "X", "ntrials", "nsubs", "gender")
  params <- c("mu_a_rew", "alpha_a_rew",
              "mu_a_pun", "alpha_a_pun",
              "mu_K",
              "mu_omega_f", "alpha_omega_f",
              "mu_omega_p",
              "lambda_a_rew", "lambda_a_pun", "lambda_K",
              "lambda_omega_f", "lambda_omega_p")
  
  # Run the model in parallel (using 3 chains and clusters)
  samples <- jags.parallel(data = data_jags, inits = NULL, params,
                           model.file = "hier_ORL_gender_new.txt", 
                           n.chains = 3, n.iter = 4000, n.burnin = 1000, 
                           n.thin = 1, n.cluster = 3)
  
  # Save true parameters for this iteration
  true_mu_a_rew[i]   <- mu_a_rew
  true_mu_a_pun[i]   <- mu_a_pun
  true_mu_K[i]       <- mu_K
  true_mu_omega_f[i] <- mu_omega_f
  true_mu_omega_p[i] <- mu_omega_p
  
  true_alpha_a_rew[i]   <- alpha_a_rew
  true_alpha_a_pun[i]   <- alpha_a_pun
  true_alpha_omega_f[i] <- alpha_omega_f
  
  # Extract inferred parameters (using the maximum of the posterior density)
  Y <- samples$BUGSoutput$sims.list
  infer_mu_a_rew[i]   <- MPD(Y$mu_a_rew)
  infer_mu_a_pun[i]   <- MPD(Y$mu_a_pun)
  infer_mu_K[i]       <- MPD(Y$mu_K)
  infer_mu_omega_f[i] <- MPD(Y$mu_omega_f)
  infer_mu_omega_p[i] <- MPD(Y$mu_omega_p)
  
  true_lambda_a_rew[i]   <- 1 / (sigma_a_rew^2)
  true_lambda_a_pun[i]   <- 1 / (sigma_a_pun^2)
  true_lambda_K[i]       <- 1 / (sigma_K^2)
  true_lambda_omega_f[i] <- 1 / (sigma_omega_f^2)
  true_lambda_omega_p[i] <- 1 / (sigma_omega_p^2)
  
  infer_lambda_a_rew[i]   <- MPD(Y$lambda_a_rew)
  infer_lambda_a_pun[i]   <- MPD(Y$lambda_a_pun)
  infer_lambda_K[i]       <- MPD(Y$lambda_K)
  infer_lambda_omega_f[i] <- MPD(Y$lambda_omega_f)
  infer_lambda_omega_p[i] <- MPD(Y$lambda_omega_p)
  
  infer_alpha_a_rew[i]   <- MPD(Y$alpha_a_rew)
  infer_alpha_a_pun[i]   <- MPD(Y$alpha_a_pun)
  infer_alpha_omega_f[i] <- MPD(Y$alpha_omega_f)
  
  # -----------------------------
  # Print progress information
  # -----------------------------
  progress_pct <- round(i / niterations * 100, 2)
  cat("Progress:", progress_pct, "% complete.\n")
  elapsed <- Sys.time() - start_time
  cat("Time elapsed since start:", elapsed, "\n\n")
  
  # -----------------------------
  # Save intermediate progress to file (so that progress is not lost if interrupted)
  # -----------------------------
  progress_file <- paste0("hier_ORL_recovery_gender_latest_seed_", seed, ".RData")
  save(ntrials_all, niterations,
       true_mu_a_rew, infer_mu_a_rew,
       true_mu_a_pun, infer_mu_a_pun,
       true_mu_K, infer_mu_K,
       true_mu_omega_f, infer_mu_omega_f,
       true_mu_omega_p, infer_mu_omega_p,
       true_lambda_a_rew, infer_lambda_a_rew,
       true_lambda_a_pun, infer_lambda_a_pun,
       true_lambda_K, infer_lambda_K,
       true_lambda_omega_f, infer_lambda_omega_f,
       true_lambda_omega_p, infer_lambda_omega_p,
       true_alpha_a_rew, infer_alpha_a_rew,
       true_alpha_a_pun, infer_alpha_a_pun,
       true_alpha_omega_f, infer_alpha_omega_f,
       file = progress_file)
}
end_time <- Sys.time()
cat("Total time elapsed:", end_time - start_time, "\n")

# -----------------------------
# 7. Save the final results
# -----------------------------
save(ntrials_all, niterations,
     # Mu parameters
     true_mu_a_rew, infer_mu_a_rew,
     true_mu_a_pun, infer_mu_a_pun,
     true_mu_K, infer_mu_K,
     true_mu_omega_f, infer_mu_omega_f,
     true_mu_omega_p, infer_mu_omega_p,
     # Lambda parameters
     true_lambda_a_rew, infer_lambda_a_rew,
     true_lambda_a_pun, infer_lambda_a_pun,
     true_lambda_K, infer_lambda_K,
     true_lambda_omega_f, infer_lambda_omega_f,
     true_lambda_omega_p, infer_lambda_omega_p,
     # Gender effects (alpha)
     true_alpha_a_rew, infer_alpha_a_rew,
     true_alpha_a_pun, infer_alpha_a_pun,
     true_alpha_omega_f, infer_alpha_omega_f,
     file = "hier_ORL_recovery_gender_latest.RData")

# -----------------------------
# 8. Plot traceplots of the last model run (optional)
# -----------------------------
par(mar = c(4, 4, 2, 1))
traceplot(samples, mfrow = c(7, 2))