# Load packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# Load results file
setwd('/work/LinusJoakimBackström#7558/DM exam')
load("results/ORL_gender_latest_results.RData")

# Extract posterior samples
post <- fit$BUGSoutput$sims.list

# Posterior probabilities that females have higher 
# values for the corresponding parameter
p_alpha_a_rew <- mean(post$alpha_a_rew > 0)
p_alpha_a_pun <- mean(post$alpha_a_pun > 0)
p_alpha_omega_f <- mean(post$alpha_omega_f > 0)

# Extract alpha posteriors
alpha_data <- data.frame(
  Parameter = rep(c("alpha_a_rew", "alpha_a_pun", "alpha_omega_f"), each = length(post$alpha_a_rew)),
  Value = c(post$alpha_a_rew, post$alpha_a_pun, post$alpha_omega_f)
)

# Summary table for numerical reporting
alpha_summary <- alpha_data %>%
  group_by(Parameter) %>%
  summarise(
    Mean = mean(Value),
    CI_lower = quantile(Value, 0.025),
    CI_upper = quantile(Value, 0.975),
  ) %>%
  ungroup()


# --- POSTERIOR PLOTS ---

# Density data with CI bounds (for plotting)
plot_data <- alpha_data %>%
  group_by(Parameter) %>%
  do({
    dens <- density(.$Value, adjust = 1.5)
    df <- tibble(x = dens$x, y = dens$y)
    ci <- quantile(.$Value, c(0.025, 0.975))
    df$ci_low <- ci[1]
    df$ci_high <- ci[2]
    df
  }) %>%
  ungroup()

# Second df for the shaded ribbon
ribbon_data <- plot_data %>%
  filter(x >= ci_low & x <= ci_high)

# Plot Posterior Distributions of Gender Differences (alpha)
ggplot(plot_data, aes(x = x, y = y)) +
  
  # Draw the full density curve as a line
  geom_line(color = "black") +
  # Add shaded ribbon for 95% CI
  geom_ribbon(data = ribbon_data, aes(ymin = 0, ymax = y), 
              fill = "skyblue", alpha = 0.75) +
  # Add vertical line at zero (null hypothesis)
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = "red", linewidth = 0.5) +
  facet_wrap(~Parameter, scales = "free") +
  
  # Labels and theme
  labs(x = "Difference (Female - Male)",
       y = "Density") +
  theme_minimal()


# Extract posteriors for each group
female_draws <- data.frame(
  a_rew = post$mu_a_rew + 0.5 * post$alpha_a_rew,
  a_pun = post$mu_a_pun + 0.5 * post$alpha_a_pun,
  omega_f = post$mu_omega_f + 0.5 * post$alpha_omega_f,
  group = "Female"
)
male_draws <- data.frame(
  a_rew = post$mu_a_rew - 0.5 * post$alpha_a_rew,
  a_pun = post$mu_a_pun - 0.5 * post$alpha_a_pun,
  omega_f = post$mu_omega_f - 0.5 * post$alpha_omega_f,
  group = "Male"
)

plot_data <- bind_rows(female_draws, male_draws)
plot_long <- pivot_longer(plot_data, cols = -group, names_to = "parameter", values_to = "value")

# Compute density
density_data <- plot_long %>%
  group_by(group, parameter) %>%
  do({
    d <- density(.$value, adjust = 1)
    tibble(x = d$x, y = d$y)
  })

# Compute 95% CIs
ci_bounds <- plot_long %>%
  group_by(group, parameter) %>%
  summarise(
    ci_low = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    .groups = "drop"
  )

density_with_ci <- density_data %>%
  left_join(ci_bounds, by = c("group", "parameter")) %>%
  mutate(in_ci = x >= ci_low & x <= ci_high)

# Plot
ggplot() +
  geom_ribbon(data = density_with_ci %>% filter(in_ci),
              aes(x = x, ymin = 0, ymax = y, fill = group),
              alpha = 0.5) +
  geom_line(data = density_with_ci,
            aes(x = x, y = y, color = group)) +
  facet_wrap(~parameter, scales = "free") +
  labs(x = "Parameter Value", y = "Density") +
  theme_minimal() +
  theme(strip.text = element_text(size = 11))



