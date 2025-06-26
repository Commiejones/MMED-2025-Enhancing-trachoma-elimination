# Load required libraries
library(readr)
library(dplyr)
library(deSolve)
library(tidyr)
library(ggplot2)

# ─────────────────────────────────────────────────────────────
# Load and preprocess data

individual_data <- read_csv("trachoma_sero_transmission_analysis_indiv.csv", show_col_types = FALSE)

# Filter for WUHA site and valid serology
individual_data <- individual_data %>%
  filter(location_name == "Wag Hemra, Ethiopia (WUHA)") %>%
  mutate(sero_valid = !is.na(pgp3_pos))

# Define 1-year age bins from 1 to 9 years
age_bins <- 1:10

# Calculate seroprevalence by age group
seroprev_data <- individual_data %>%
  mutate(age_group = cut(age_years, breaks = age_bins, include.lowest = TRUE, right = FALSE)) %>%
  filter(sero_valid) %>%
  group_by(age_group) %>%
  summarise(
    n_sero = n(),
    seropositive = sum(pgp3_pos == 1),
    seroprev_obs = seropositive / n_sero,
    .groups = "drop"
  ) %>%
  mutate(
    age_mid = (as.numeric(sub("\\[", "", sapply(strsplit(as.character(age_group), ","), `[`, 1))) +
                 as.numeric(sub("\\)", "", sapply(strsplit(as.character(age_group), ","), `[`, 2)))) / 2,
    var_seroprev = seroprev_obs * (1 - seroprev_obs) / n_sero
  ) %>%
  filter(!is.na(age_mid)) %>%
  select(age_mid, n_sero, seroprev_obs, var_seroprev)

model_data <- seroprev_data

# ─────────────────────────────────────────────────────────────
# SIDR Model Functions with correct ODEs including reinfection

# Constant lambda model
sidr_model_constant <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambda <- lambda_0
    
    dS <- -lambda * S
    dI <- lambda * S - gamma * I + alpha * lambda * R
    dD <- gamma * I - rho * D
    dR <- rho * D - alpha * lambda * R
    
    list(c(dS, dI, dD, dR))
  })
}

predict_seroprev_constant <- function(lambda_0, times, parameters_fixed) {
  if (lambda_0 < 0) return(rep(NA, length(times)))
  parameters <- c(lambda_0 = lambda_0, parameters_fixed)
  state <- c(S = 1000, I = 0, D = 0, R = 0)  # all susceptible initially
  out <- tryCatch(
    ode(y = state, times = times, func = sidr_model_constant, parms = parameters),
    error = function(e) return(NULL)
  )
  if (is.null(out)) return(rep(NA, length(times)))
  out_df <- as.data.frame(out)
  seroprev <- (out_df$D + out_df$R) / rowSums(out_df[, c("S", "I", "D", "R")])
  return(seroprev)
}

loss_function_constant <- function(lambda_0, model_data, parameters_fixed) {
  pred <- predict_seroprev_constant(lambda_0, model_data$age_mid, parameters_fixed)
  if (any(is.na(pred))) return(1e10)
  residuals <- (model_data$seroprev_obs - pred)^2 / model_data$var_seroprev
  sum(residuals, na.rm = TRUE)
}

# Decreasing lambda model with exponential decay
sidr_model_decreasing <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambda <- lambda_0 * exp(-time / F0)
    
    dS <- -lambda * S
    dI <- lambda * S - gamma * I + alpha * lambda * R
    dD <- gamma * I - rho * D
    dR <- rho * D - alpha * lambda * R
    
    list(c(dS, dI, dD, dR))
  })
}

predict_seroprev_decreasing <- function(lambda_0, F0, times, parameters_fixed) {
  if (lambda_0 < 0 || F0 <= 0) return(rep(NA, length(times)))
  parameters <- c(lambda_0 = lambda_0, F0 = F0, parameters_fixed)
  state <- c(S = 1000, I = 0, D = 0, R = 0)
  out <- tryCatch(
    ode(y = state, times = times, func = sidr_model_decreasing, parms = parameters),
    error = function(e) return(NULL)
  )
  if (is.null(out)) return(rep(NA, length(times)))
  out_df <- as.data.frame(out)
  seroprev <- (out_df$D + out_df$R) / rowSums(out_df[, c("S", "I", "D", "R")])
  return(seroprev)
}

loss_function_decreasing <- function(par, model_data, parameters_fixed) {
  lambda_0 <- par[1]
  F0 <- par[2]
  pred <- predict_seroprev_decreasing(lambda_0, F0, model_data$age_mid, parameters_fixed)
  if (any(is.na(pred))) return(1e10)
  residuals <- (model_data$seroprev_obs - pred)^2 / model_data$var_seroprev
  sum(residuals, na.rm = TRUE)
}

# ─────────────────────────────────────────────────────────────
# Fit models

parameters_fixed <- c(gamma = 2, rho = 1, alpha = 0.1)  # example fixed rates per year

fit_const <- optim(
  par = 0.2,
  fn = loss_function_constant,
  model_data = model_data,
  parameters_fixed = parameters_fixed,
  method = "L-BFGS-B",
  lower = 0.001,
  upper = 5
)

fit_decr <- optim(
  par = c(0.2, 3),
  fn = loss_function_decreasing,
  model_data = model_data,
  parameters_fixed = parameters_fixed,
  method = "L-BFGS-B",
  lower = c(0.001, 0.1),
  upper = c(5, 20)
)

# ─────────────────────────────────────────────────────────────
# Generate predictions with fitted parameters

model_data <- model_data %>%
  mutate(
    predicted_constant = predict_seroprev_constant(fit_const$par, age_mid, parameters_fixed),
    predicted_decreasing = predict_seroprev_decreasing(fit_decr$par[1], fit_decr$par[2], age_mid, parameters_fixed)
  )

# ─────────────────────────────────────────────────────────────
# Plot results with legend

plot_data <- model_data %>%
  pivot_longer(
    cols = starts_with("predicted"),
    names_to = "model",
    values_to = "seroprev"
  ) %>%
  mutate(
    model = recode(model,
                   predicted_constant = "Constant λ",
                   predicted_decreasing = "Decreasing λ")
  )

ggplot() +
  geom_point(data = model_data, aes(x = age_mid, y = seroprev_obs), color = "black", size = 3, shape = 16) +
  geom_line(data = plot_data, aes(x = age_mid, y = seroprev, color = model), size = 1.3) +
  scale_color_manual(values = c("Constant λ" = "blue", "Decreasing λ" = "red")) +
  labs(
    x = "Age (years)",
    y = "Seroprevalence",
    title = "Fit of SIDR Models to Seroprevalence Data in WUHA",
    color = "Model"
  ) +
  theme_minimal(base_size = 14)

# ─────────────────────────────────────────────────────────────
# Print fitted parameters

cat("Fitted constant λ:", fit_const$par, "\n")
cat("Fitted decreasing λ parameters: lambda_0 =", fit_decr$par[1], ", F =", fit_decr$par[2], "\n")

