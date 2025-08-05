if (!require(pacman)) install.packages("pacman")
pacman::p_load(R2jags, tidyverse, ggplot2, cowplot, truncnorm)

setwd('/work/LinusJoakimBackström#7558/DM exam')
load("results/ORL_gender_lambdas_results.RData") 

gamtest = read.csv("data/Gamtest_PreprocessedOSF.csv", row.names = 1)
IGT = read.csv("data/Preprocessed_IGT_OSF.csv", row.names = 1)

# Preprocess data exactly as in the model fitting script
gamtest <- gamtest[gamtest$Gender != "Other", ]
IGT <- IGT[IGT$ID %in% gamtest$ID, ]
gamtest$GenderBinary <- ifelse(gamtest$Gender == "Female", 1, 0)
IGT <- merge(IGT, gamtest[, c("ID", "GenderBinary")], by = "ID")

choice_matrix <- IGT %>%
  select(ID, Trial, choice) %>%
  pivot_wider(names_from = ID, values_from = choice) %>%
  select(-Trial) %>%
  as.matrix()

outcome_matrix <- IGT %>%
  select(ID, Trial, net_outcome) %>%
  mutate(net_outcome = net_outcome / 100) %>% # scale outcomes
  pivot_wider(names_from = ID, values_from = net_outcome) %>%
  select(-Trial) %>%
  as.matrix()

gender_vector <- gamtest %>%  
  arrange(match(ID, colnames(choice_matrix))) %>%
  pull(GenderBinary)

jags_data <- list(
  x = t(choice_matrix),
  X = t(outcome_matrix),
  nsubs = ncol(choice_matrix),
  ntrials = rep(100, ncol(choice_matrix)),
  gender = gender_vector
)

hier_ORL_gender_sim <- function(payoff, nsubs, ntrials, gender_vec,
                                mu_a_rew, alpha_a_rew,
                                mu_a_pun, alpha_a_pun,
                                mu_K,
                                mu_omega_f, alpha_omega_f,
                                mu_omega_p,
                                sigma_a_rew, sigma_a_pun, sigma_K,
                                sigma_omega_f, sigma_omega_p) {
  
  x <- array(NA, c(nsubs, ntrials))      # Choices
  X <- array(NA, c(nsubs, ntrials))      # Outcomes
  Ev <- array(NA, c(nsubs, ntrials, 4))  # Expected values
  
  for (s in 1:nsubs) {
    # Draw individual-level parameters based on the group-level posteriors
    a_rew <- rtruncnorm(1, a=0, b=1, mean = mu_a_rew + alpha_a_rew * (gender_vec[s] - 0.5), sd = sigma_a_rew)
    a_pun <- rtruncnorm(1, a=0, b=1, mean = mu_a_pun + alpha_a_pun * (gender_vec[s] - 0.5), sd = sigma_a_pun)
    K <- rtruncnorm(1, a=0, b=5, mean = mu_K, sd = sigma_K)
    omega_f <- rnorm(1, mean = mu_omega_f + alpha_omega_f * (gender_vec[s] - 0.5), sd = sigma_omega_f)
    omega_p <- rnorm(1, mean = mu_omega_p, sd = sigma_omega_p)
    theta <- 1 # Fixed
    
    # Initialize subject-specific arrays
    Ef <- matrix(0, ntrials, 4)
    PS <- matrix(1, ntrials, 4)
    
    # Initialize values for the first trial
    Ev[s, 1, ] <- 0
    
    for (t in 1:(ntrials - 1)) {
      # Calculate V and choice probabilities for trial t
      V <- Ev[s, t, ] + Ef[t, ] * omega_f + PS[t, ] * omega_p
      exp_p <- exp(theta * V)
      p <- exp_p / sum(exp_p)
      
      # Ensure probabilities are valid to prevent crashes
      p[is.na(p) | is.infinite(p)] <- 0
      if (sum(p) == 0) p <- rep(0.25, 4)
      
      # Make choice for trial t
      x[s, t] <- sample(1:4, 1, prob = p)
      X[s, t] <- payoff[t, x[s, t]]
      
      # --- Update for next trial (t+1) ---
      signX <- ifelse(X[s, t] < 0, -1, 1)
      
      for (d in 1:4) {
        # Chosen deck updates
        if (d == x[s, t]) {
          # Value Update
          Ev_update <- if (X[s, t] >= 0) {
            Ev[s, t, d] + a_rew * (X[s, t] - Ev[s, t, d])
          } else {
            Ev[s, t, d] + a_pun * (X[s, t] - Ev[s, t, d])
          }
          Ev[s, t + 1, d] <- Ev_update
          
          # Frequency Update
          Ef_cho <- if (X[s, t] >= 0) {
            Ef[t, d] + a_rew * (signX - Ef[t, d])
          } else {
            Ef[t, d] + a_pun * (signX - Ef[t, d])
          }
          Ef[t + 1, d] <- Ef_cho
          
          # Perseverance Update
          PS[t + 1, d] <- 1 / (1 + K)
          
        } else { # Unchosen deck updates
          Ev[s, t + 1, d] <- Ev[s, t, d]
          
          # Frequency Update
          Ef_not <- if (X[s, t] >= 0) {
            Ef[t, d] + a_pun * ((-signX / 3) - Ef[t, d])
          } else {
            Ef[t, d] + a_rew * ((-signX / 3) - Ef[t, d])
          }
          Ef[t + 1, d] <- Ef_not
          
          # Perseverance Update
          PS[t + 1, d] <- PS[t, d] / (1 + K)
        }
      }
    }
    # Make final choice on the last trial to avoid index out of bounds error
    V_last <- Ev[s, ntrials, ] + Ef[ntrials, ] * omega_f + PS[ntrials, ] * omega_p
    exp_p_last <- exp(theta * V_last)
    p_last <- exp_p_last / sum(exp_p_last)
    p_last[is.na(p_last) | is.infinite(p_last)] <- 0
    if (sum(p_last) == 0) p_last <- rep(0.25, 4)
    x[s, ntrials] <- sample(1:4, 1, prob = p_last)
  }
  
  return(list(x = x))
}

post <- fit$BUGSoutput$sims.list
n_iter <- fit$BUGSoutput$n.sims
n_subs <- jags_data$nsubs
n_trials <- jags_data$ntrials[1]
gender_vector <- jags_data$gender

# Load the IGT payoff structure
payoff_matrix <- readRDS("payoff_structure.rds")

n_sim <- 500 # Number of posterior samples to simulate. 500 is a good number.
sim_choices <- array(NA, c(n_sim, n_subs, n_trials))

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
for (i in 1:n_sim) {
  # Randomly pick one iteration from the posterior
  iter_idx <- sample(1:n_iter, 1)
  
  sigma_a_rew <- 1 / sqrt(post$lambda_a_rew[iter_idx])
  sigma_a_pun <- 1 / sqrt(post$lambda_a_pun[iter_idx])
  sigma_K <- 1 / sqrt(post$lambda_K[iter_idx])
  sigma_omega_f <- 1 / sqrt(post$lambda_omega_f[iter_idx])
  sigma_omega_p <- 1 / sqrt(post$lambda_omega_p[iter_idx])
  
  sim_data <- hier_ORL_gender_sim(
    payoff = payoff_matrix,
    nsubs = n_subs,
    ntrials = n_trials,
    gender_vec = gender_vector,
    mu_a_rew = post$mu_a_rew[iter_idx], alpha_a_rew = post$alpha_a_rew[iter_idx],
    mu_a_pun = post$mu_a_pun[iter_idx], alpha_a_pun = post$alpha_a_pun[iter_idx],
    mu_K = post$mu_K[iter_idx],
    mu_omega_f = post$mu_omega_f[iter_idx], alpha_omega_f = post$alpha_omega_f[iter_idx],
    mu_omega_p = post$mu_omega_p[iter_idx],
    sigma_a_rew = sigma_a_rew, sigma_a_pun = sigma_a_pun, sigma_K = sigma_K,
    sigma_omega_f = sigma_omega_f, sigma_omega_p = sigma_omega_p
  )
  
  sim_choices[i, , ] <- sim_data$x
  setTxtProgressBar(pb, i)
}
close(pb)


# Function to calculate net score over blocks
calculate_net_score <- function(choice_data, n_blocks = 5) {
  n_trials <- ncol(choice_data)
  block_size <- n_trials / n_blocks
  
  # Decks 1 & 2 (A, B) are disadvantageous, 3 & 4 (C, D) are advantageous
  adv_choices <- (choice_data == 3) | (choice_data == 4)
  disadv_choices <- (choice_data == 1) | (choice_data == 2)
  
  net_score <- adv_choices - disadv_choices
  
  # Calculate mean net score per block
  block_scores <- apply(net_score, 1, function(subj_scores) {
    sapply(1:n_blocks, function(b) {
      start <- (b - 1) * block_size + 1
      end <- b * block_size
      mean(subj_scores[start:end], na.rm = TRUE)
    })
  })
  
  return(t(block_scores)) # Return as subjects x blocks
}

# Calculate net score for REAL data
real_net_scores <- calculate_net_score(t(choice_matrix))
real_mean_scores <- apply(real_net_scores, 2, mean)
real_sem_scores <- apply(real_net_scores, 2, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))

# Calculate net score for SIMULATED data
sim_mean_scores <- array(NA, c(n_sim, 5))
for (i in 1:n_sim) {
  sim_net_scores <- calculate_net_score(sim_choices[i, , ])
  sim_mean_scores[i, ] <- apply(sim_net_scores, 2, mean)
}

# Get the mean and 95% CI for the simulated scores
sim_final_mean <- apply(sim_mean_scores, 2, mean)
sim_final_ci <- apply(sim_mean_scores, 2, function(x) quantile(x, c(0.025, 0.975)))

# Create data frames for plotting with ggplot
plot_df_real <- data.frame(
  block = 1:5,
  mean_score = real_mean_scores,
  upper = real_mean_scores + real_sem_scores,
  lower = real_mean_scores - real_sem_scores,
  type = "Observed Data"
)

plot_df_sim <- data.frame(
  block = 1:5,
  mean_score = sim_final_mean,
  upper = sim_final_ci[2,],
  lower = sim_final_ci[1,],
  type = "Model Prediction"
)

# PLOT
# Posterior Predictive Check: Model vs. Real Data
# Net Score (Advantageous - Disadvantageous Choices) Across Trial Blocks
ppc_plot <- ggplot() +
  geom_ribbon(data = plot_df_sim, aes(x = block, ymin = lower, ymax = upper, fill = type), alpha = 0.5) +
  geom_line(data = plot_df_sim, aes(x = block, y = mean_score, color = type), linetype = "dashed", linewidth = 1) +
  geom_line(data = plot_df_real, aes(x = block, y = mean_score, color = type), linewidth = 1) +
  geom_point(data = plot_df_real, aes(x = block, y = mean_score, color = type), size = 3) +
  geom_errorbar(data = plot_df_real, aes(x = block, ymin = lower, ymax = upper), width = 0.1, linewidth = 0.7, color = "black") +
  scale_color_manual(name = "Legend", values = c("Model Prediction" = "deepskyblue4", "Observed Data" = "black")) +
  scale_fill_manual(name = "Legend", values = c("Model Prediction" = "skyblue")) +
  scale_x_continuous(breaks = 1:5) +
  labs(x = "Trial Block", y = "Mean Net Score") +
  theme_minimal() +
  guides(fill = "none") +
  theme(legend.position = "top") +
  theme(legend.title=element_blank())

print(ppc_plot)