install.packages("pacman")
pacman::p_load(ggplot2, tidyverse)

# Plot 1: Recovery for Gender Effects (alpha)
recovery_long_alpha <- recovery_data %>%
  select(starts_with("true_alpha"), starts_with("infer_alpha")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "parameter"),
    names_pattern = "(true|infer)_(alpha_.*)"
  ) %>%
  mutate(parameter_label = case_when(
    parameter == "alpha_a_rew" ~ "alpha[A[rew]]",
    parameter == "alpha_a_pun" ~ "alpha[A[pun]]",
    TRUE ~ "alpha[omega[f]]"
  ))

ggplot(recovery_long_alpha, aes(x = true, y = infer)) +
  geom_point(alpha = 0.4, color = "dodgerblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~parameter_label, scales = "free", labeller = label_parsed) +
  labs(x = "True Value", y = "Inferred Value") +
  theme_bw() +
  theme(strip.text = element_text(size = 14))

# Plot 2: Recovery for Group Means (mu)
recovery_long_mu <- recovery_data %>%
  select(starts_with("true_mu"), starts_with("infer_mu")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "parameter"),
    names_pattern = "(true|infer)_(mu_.*)"
  ) %>%
  mutate(parameter_label = case_when(
    parameter == "mu_a_rew"   ~ "mu[A[rew]]",
    parameter == "mu_a_pun"   ~ "mu[A[pun]]",
    parameter == "mu_K"       ~ "mu[K]",
    parameter == "mu_omega_f" ~ "mu[omega[f]]",
    TRUE                     ~ "mu[omega[p]]"
  ))

ggplot(recovery_long_mu, aes(x = true, y = infer)) +
  geom_point(alpha = 0.4, color = "dodgerblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~parameter_label, scales = "free", labeller = label_parsed) +
  labs(x = "True Value", y = "Inferred Value") +
  theme_bw() +
  theme(strip.text = element_text(size = 14))