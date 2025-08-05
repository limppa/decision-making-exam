# Load packages
install.packages("pacman")
pacman::p_load(R2jags, tidyverse)

# Set working directory
setwd('/work/LinusJoakimBackström#7558/DM exam')

# Load data
gamtest = read.csv("data/Gamtest_PreprocessedOSF.csv", row.names = 1)
IGT = read.csv("data/Preprocessed_IGT_OSF.csv", row.names = 1)

# Exclude one participant with "Other" gender
gamtest <- gamtest[gamtest$Gender != "Other", ]
IGT <- IGT[IGT$ID %in% gamtest$ID, ]

# Create binary gender column for gender model (1=Female, 0=Male)
gamtest$GenderBinary <- ifelse(gamtest$Gender == "Female", 1, 0)
# Add gender info to IGT df
IGT <- merge(IGT, gamtest[, c("ID", "GenderBinary")], by = "ID")

# Create matrices
choice_matrix <- IGT %>%
  select(ID, Trial, choice) %>%
  pivot_wider(names_from = ID, values_from = choice) %>%
  select(-Trial) %>%
  as.matrix()

outcome_matrix <- IGT %>%
  select(ID, Trial, net_outcome) %>%
  mutate(net_outcome = net_outcome / 100) %>% # scale outcomes to match
  pivot_wider(names_from = ID, values_from = net_outcome) %>%
  select(-Trial) %>%
  as.matrix()

# Gender vector
gender_vector <- gamtest %>% 
  arrange(match(ID, colnames(choice_matrix))) %>%
  pull(GenderBinary)

# JAGS data list
jags_data <- list(
  x = t(choice_matrix),
  X = t(outcome_matrix),
  nsubs = ncol(choice_matrix),
  ntrials = rep(100, ncol(choice_matrix)),
  gender = gender_vector
)

# Parameters to monitor
params <- c("mu_a_rew", "alpha_a_rew",
            "mu_a_pun", "alpha_a_pun",
            "mu_omega_f", "alpha_omega_f",
            "mu_K", "mu_omega_p")

# Run model
fit <- jags.parallel(
  data = jags_data,
  inits = NULL,
  parameters.to.save = params,
  model.file = "hier_ORL_gender_latest.txt",
  n.chains = 4,
  n.iter = 10000,
  n.burnin = 3000,  
  n.thin = 2,
  n.cluster = 4
)

# Save results
save(fit, file = "results/ORL_gender_latest_results.RData")