library(readr)
library(dplyr)
library(deSolve)
library(tidyr)
library(ggplot2)

# ─────────────────────────────────────────────────────────────
# Load and preprocess data (same as before)

individual_data <- read_csv("trachoma_sero_transmission_analysis_indiv.csv", show_col_types = FALSE)

individual_data <- individual_data %>%
  filter(location_name == "Wag Hemra, Ethiopia (WUHA)", year == 2016) %>%
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
# Corrected SIDR model with logistic lambda(t)

sidr_model_logistic_corrected <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Logistic force of infection
    lambda <- lambda_0 / (1 + exp((time - t_half) / k))
    
    dS <- -lambda * S
    dI <- lambda * S - gamma * I + alpha * lambda * R
    dD <- gamma * I - rho * D
    dR <- rho * D - alpha * lambda * R
    
    list(c(dS, dI, dD, dR))
  })
}

predict_seroprev_logistic_corrected <- function(lambda_0, t_half, k, alpha, times, parameters_fixed) {
  if (lambda_0 < 0 || k <= 0 || alpha < 0) return(rep(NA, length(times)))
  parameters <- c(lambda_0 = lambda_0, t_half = t_half, k = k, alpha = alpha, parameters_fixed)
  state <- c(S = 1000, I = 0, D = 0, R = 0)
  
  out <- tryCatch(
    ode(y = state, times = times, func = sidr_model_logistic_corrected, parms = parameters),
    error = function(e) return(NULL)
  )
  
  if (is.null(out)) return(rep(NA, length(times)))
  out_df <- as.data.frame(out)
  seroprev <- (out_df$D + out_df$R) / rowSums(out_df[, c("S", "I", "D", "R")])
  return(seroprev)
}

loss_function_logistic_corrected <- function(par, model_data, parameters_fixed) {
  lambda_0 <- par[1]
  t_half <- par[2]
  k <- par[3]
  alpha <- par[4]
  
  pred <- predict_seroprev_logistic_corrected(lambda_0, t_half, k, alpha, model_data$age_mid, parameters_fixed)
  if (any(is.na(pred))) return(1e10)
  
  residuals <- (model_data$seroprev_obs - pred)^2 / model_data$var_seroprev
  sum(residuals, na.rm = TRUE)
}

# ─────────────────────────────────────────────────────────────
# Fit the corrected model

parameters_fixed <- c(gamma = 2, rho = 1, sigma = 0.5)  # fixed rates per year

fit_logistic_corrected <- optim(
  par = c(0.3, 4, 1, 0.1),  # initial guesses: lambda_0, t_half, k, alpha
  fn = loss_function_logistic_corrected,
  model_data = model_data,
  parameters_fixed = parameters_fixed,
  method = "L-BFGS-B",
  lower = c(0.001, 0, 0.1, 0),
  upper = c(5, 10, 10, 5)
)

# ─────────────────────────────────────────────────────────────
# Generate predictions with fitted parameters

model_data <- model_data %>%
  mutate(
    predicted_logistic_corrected = predict_seroprev_logistic_corrected(
      fit_logistic_corrected$par[1],
      fit_logistic_corrected$par[2],
      fit_logistic_corrected$par[3],
      fit_logistic_corrected$par[4],
      age_mid,
      parameters_fixed
    )
  )

# ─────────────────────────────────────────────────────────────
# Plot results

ggplot() +
  geom_point(data = model_data, aes(x = age_mid, y = seroprev_obs), color = "black", size = 3) +
  geom_line(data = model_data, aes(x = age_mid, y = predicted_logistic_corrected), color = "darkgreen", size = 1.3) +
  labs(
    x = "Age (years)",
    y = "Seroprevalence",
    title = "Fit of Corrected SIDR Model with Logistic λ(t) to Seroprevalence Data",
    color = "Model"
  ) +
  theme_minimal(base_size = 14)

# ─────────────────────────────────────────────────────────────
# Print fitted parameters

cat("Fitted parameters:\n")
cat("lambda_0 =", fit_logistic_corrected$par[1], "\n")
cat("t_half =", fit_logistic_corrected$par[2], "\n")
cat("k =", fit_logistic_corrected$par[3], "\n")
cat("alpha =", fit_logistic_corrected$par[4], "\n")




