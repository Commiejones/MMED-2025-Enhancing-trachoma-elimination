################################################################################
## PACKAGES AND LIBRARIES ##
################################################################################
library(readr)       
library(tidyverse)   
library(dplyr)       
library(ggplot2)     
library(purrr)       
library(broom)       

################################################################################
## LOAD DATA ##
################################################################################
cluster_data <- read_csv("Documents/AIMS/MMED/MMED R/MMED TRACHOMA/trachoma_sero_transmission_analysis_cluster.csv")
indiv_data <- read_csv("Documents/AIMS/MMED/MMED R/MMED TRACHOMA/trachoma_sero_transmission_analysis_indiv.csv")
study_data <- read_csv("Documents/AIMS/MMED/MMED R/MMED TRACHOMA/trachoma_sero_transmission_analysis_study.csv")

################################################################################
## BASIC ANALYSIS ##
################################################################################
summary(cluster_data)   # Summarize cluster-level data
summary(indiv_data)     # Summarize individual-level data
summary(study_data)     # Summarize study-level data

# str() or head()

################################################################################
## PREPARE AND ANALYZE INDIVIDUAL-LEVEL SEROPREVALENCE ##
################################################################################
unique(indiv_data$age_years)  #Ages are only from 1-9

indiv <- indiv_data %>% 
  filter(!is.na(pgp3_pos)) %>%                   # Exclude rows with missing serology results
  mutate(age = as.numeric(age_years)) %>%        # Ensure age is numeric
  filter(age >= 1 & age <= 9) %>%                # Keep ages 1 to 9
  group_by(age) %>%                              # Group by age
  summarise(n = n(),                             # Total number of individuals
            pos = sum(pgp3_pos),                 # Number seropositive
            seroprev = pos / n)                  # Calculate seroprevalence

# Fit serocatalytic model (no seroreversion)
fit <- nls(seroprev ~ 1 - exp(-lambda * age), #Non-linear model fitting
           data = indiv,
           start = list(lambda = 0.1))

summary(fit)    
# View model summary
confint(fit)
lambda_est <- coef(fit)[["lambda"]]               # Extract lambda estimate

# Plot fitted model on normal scale
ggplot(indiv, aes(x = age, y = seroprev)) +
  geom_point(size = 3) +
  stat_function(fun = function(a) 1 - exp(-lambda_est * a),
                color = "blue", size = 1) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  labs(title = "Fitted Seroprevalence Model",
       y = "Seroprevalence", x = "Age")

# Plot model on log-odds scale
indiv$seroprev_clipped <- pmin(pmax(indiv$seroprev, 0.001), 0.999)

ggplot(indiv, aes(x = age, y = seroprev_clipped)) +
  geom_point(size = 3) +
  stat_function(fun = function(a) lambda_est * a,  # Logit scale: logit(p) = lambda * age
                color = "blue", size = 1) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(trans = "logit",
                     breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Fitted Seroprevalence Model (Log-Odds Scale)",
       y = "Seroprevalence (Log-Odds)", x = "Age")

################################################################################
## FIT MODEL BY LOCATION ##
################################################################################

# Clean and filter individual-level data
indiv_clean <- indiv_data %>%
  filter(!is.na(pgp3_pos)) %>%              # Remove rows with missing pgp3 serology
  mutate(age = as.numeric(age_years)) %>%   # Convert age to numeric (in case it is a factor/character)
  filter(age >= 1 & age <= 9)               # Focus on children aged 1 to 9 years

# Nest data by location for group-wise modeling
cluster_nest <- indiv_clean %>%
  group_by(location_name) %>%  # Group data by location name
  nest()                       # Create a nested tibble where each row contains a sub-dataframe for one location

# Define a function to fit the serocatalytic model within each location
fit_lambda <- function(df) {
  age_summary <- df %>%
    group_by(age) %>%                              # Group by age
    summarise(n = n(),                             # Count of individuals
              pos = sum(pgp3_pos),                 # Sum of seropositives
              seroprev = pos / n,                  # Calculate seroprevalence
              .groups = "drop")
  
  # Fit the reversible catalytic model using nonlinear least squares
  tryCatch({
    model <- nls(seroprev ~ 1 - exp(-lambda * age),  # Serocatalytic model: seroprev = 1 - exp(-λ * age)
                 data = age_summary,
                 start = list(lambda = 0.1))         # Starting value for λ
    tibble(lambda = coef(model)[["lambda"]])         # Return λ as a tibble
  }, error = function(e) {
    tibble(lambda = NA)                              # Return NA if model fails to converge
  })
}

# Apply model to each nested location dataset
lambda_per_cluster <- cluster_nest %>%
  mutate(fit = map(data, fit_lambda)) %>%  # Apply `fit_lambda` to each nested data frame
  unnest(fit)                              # Unnest to get a flat tibble with one row per location

# Drop nested data column (only keep location and lambda)
lambda_per_location <- lambda_per_cluster %>%
  select(-data)

# Plot estimated λ for each location
lambda_per_location %>%
  arrange(desc(lambda)) %>%                             # Sort by descending λ
  mutate(location_name = factor(location_name, levels = location_name)) %>%  # Keep original order for plotting
  ggplot(aes(x = location_name, y = lambda)) +          # Set x and y aesthetics
  geom_col(fill = "blue") +                        # Create bar chart
  geom_text(aes(label = round(lambda, 3)),              # Add λ values as text above bars
            vjust = -0.5, size = 3.5) +
  labs(title = "Estimated Force of Infection (\u03bb) by Location",  # Title with Greek lambda symbol
       x = "Location", y = "Estimated \u03bb") +
  theme_minimal(base_size = 12) +                       # Clean minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability

################################################################################
## BOOTSTRAP CONFIDENCE INTERVALS FOR LAMBDA ##
################################################################################

# Define a bootstrap function to compute CI around λ
boot_lambda <- function(df, n = 100) {
  boot_vals <- replicate(n, {
    sample_df <- df[sample(1:nrow(df), replace = TRUE), ]  # Resample data with replacement
    fit_lambda(sample_df)$lambda                            # Get lambda from resampled data
  })
  tibble(
    lambda = mean(boot_vals, na.rm = TRUE),                # Mean λ
    lower = quantile(boot_vals, 0.025, na.rm = TRUE),      # Lower 2.5% CI
    upper = quantile(boot_vals, 0.975, na.rm = TRUE)       # Upper 97.5% CI
  )
}

# Compute bootstrap estimates for each location
lambda_ci <- indiv_clean %>%
  group_by(location_name) %>%             # Group by location
  group_modify(~ boot_lambda(.x, n = 200)) %>%  # Apply boot_lambda to each location
  ungroup()

# Plot λ with 95% confidence intervals
ggplot(lambda_ci, aes(x = reorder(location_name, lambda), y = lambda)) +
  geom_col(fill = "blue", width = 0.7) +                     # Bar plot of λ
  geom_errorbar(aes(ymin = lower, ymax = upper),                # Add error bars for CIs
                width = 0.2, color = "red") +
  geom_text(aes(label = round(lambda, 2)),                      # Add λ as text above bars
            vjust = -0.5, size = 3) +
  labs(title = "Estimated Force of Infection (\u03bb) by Location",  # Title
       x = "Location", y = "\u03bb (per year)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # Rotate and size x-axis labels
    axis.text.y = element_text(size = 10)                          # Set size for y-axis labels
  )
