# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(deSolve)

# ─────────────────────────────────────────────────────────────
# Load and preprocess data
individual_data <- read_csv("trachoma_sero_transmission_analysis_indiv.csv", show_col_types = FALSE)

individual_data <- individual_data %>%
  filter(location_name == "Wag Hemra, Ethiopia (WUHA)") %>%
  mutate(sero_valid = !is.na(pgp3_pos))

age_bins <- 1:10

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
# SIDR Model with correct equations and constant lambda

sidr_model_constant <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambda <- lambda_0  # constant force of infection
    
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
  state <- c(S = 1000, I = 0, D = 0, R = 0)  # start fully susceptible
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

# ─────────────────────────────────────────────────────────────
# SIDR Model with correct equations and decreasing lambda

sidr_model_decreasing <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    tau <- (1 / lambda_0) + (1 / gamma)
    lambda <- lambda_0 * exp(-time / tau)
    
    dS <- -lambda * S
    dI <- lambda * S - gamma * I + alpha * lambda * R
    dD <- gamma * I - rho * D
    dR <- rho * D - alpha * lambda * R
    
    list(c(dS, dI, dD, dR))
  })
}

predict_seroprev_decreasing <- function(lambda_0, times, parameters_fixed) {
  if (lambda_0 < 0) return(rep(NA, length(times)))
  parameters <- c(lambda_0 = lambda_0, parameters_fixed)
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

loss_function_decreasing <- function(lambda_0, model_data, parameters_fixed) {
  pred <- predict_seroprev_decreasing(lambda_0, model_data$age_mid, parameters_fixed)
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
  par = 0.2,
  fn = loss_function_decreasing,
  model_data = model_data,
  parameters_fixed = parameters_fixed,
  method = "L-BFGS-B",
  lower = 0.001,
  upper = 5
)

# ─────────────────────────────────────────────────────────────
# Generate predictions

model_data <- model_data %>%
  mutate(
    predicted_constant = predict_seroprev_constant(fit_const$par, age_mid, parameters_fixed),
    predicted_decreasing = predict_seroprev_decreasing(fit_decr$par, age_mid, parameters_fixed)
  )

# ─────────────────────────────────────────────────────────────
# Plot results with legend

library(tidyr)

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
    title = "SIDR Model fit to Seroprevalence Data in WUHA, Ethiopia(2016-2019)",
    subtitle = "Comparison of Constant and Decreasing Force of Infection Models",
    color = "Model"
  ) +
  theme_minimal(base_size = 14)



###
# Number of data points
n <- nrow(model_data)

# For constant lambda (1 parameter)
k_const <- 1
rss_const <- loss_function_constant(fit_const$par, model_data, parameters_fixed)
aic_const <- n * log(rss_const / n) + 2 * k_const

# For decreasing lambda (2 parameters)
k_decr <- 2
rss_decr <- loss_function_decreasing(fit_decr$par, model_data, parameters_fixed)
aic_decr <- n * log(rss_decr / n) + 2 * k_decr

cat("AIC constant λ:", aic_const, "\n")
cat("AIC decreasing λ:", aic_decr, "\n")


