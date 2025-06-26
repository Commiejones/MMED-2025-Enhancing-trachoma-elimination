################################################################################
## PACKAGES AND LIBRARIES ##
################################################################################
library(readr)       
library(dplyr)       
library(tidyr)       
library(purrr)       
library(ggplot2)
library(stringr)

################################################################################
## LOAD DATA ##
################################################################################
cluster_data <- read_csv("trachoma_sero_transmission_analysis_cluster.csv")
indiv_data <- read_csv("trachoma_sero_transmission_analysis_indiv.csv")
study_data <- read_csv("trachoma_sero_transmission_analysis_study.csv")
################################################################################
## NEW DATA ##
################################################################################
cluster_selected <- cluster_data %>%
  select(cluster_id, pcr_n_pos, pcr_n_tested, pcr_prev, tf_n_pos, tf_n_tested, 
         tf_prev, ti_n_pos, ti_n_tested, ti_prev, pgp3_n_pos, pgp3_n_tested, pgp3_prev)

indiv_selected <- indiv_data %>%
  select(year, age_years,location_name, cluster_id)

trachoma <- indiv_selected %>%
  left_join(cluster_selected, by = "cluster_id")

summary(trachoma)

################################################################################
## PLOTS ##
################################################################################

age_bins <- trachoma %>%
  mutate(age_years = case_when(
    age_years >= 1 & age_years <= 4 ~ "1-4",
    age_years >= 5 & age_years <= 9 ~ "5-9",
    TRUE ~ NA_character_  #in case there are ages outside 1-9
  ))

summary(age_bins)


## PCR POSITIVE
trachoma_pcr <- age_bins %>% 
  select(year, age_years, location_name, pcr_n_pos, pcr_n_tested, pcr_prev) %>%
  mutate(
    country = str_extract(location_name, ",\\s*[^\\(]+") %>%
      str_replace_all("^,\\s*", "")
  )
unique(trachoma_pcr$country)

pcr_by_age_country <- trachoma_pcr %>%
  group_by(year, country, age_years) %>%
  summarise(
    mean_pcr_prev = mean(pcr_prev, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(pcr_by_age_country, aes(x = age_years, y = mean_pcr_prev, color = as.factor(year))) +
  geom_line() +
  geom_point() +
  facet_wrap(~ country) +
  labs(
    title = "PCR Prevalence by Age and Year",
    x = "Age (years)",
    y = "Mean PCR Prevalence (%)",
    color = "Year"
  ) +
  theme_minimal()

#Pick year? Country? Add % values to points? Mean PCR Prevalence: rounding? NaN = 0?


## TF POSITIVE
trachoma_tf <- age_bins %>% 
  select(year, age_years, location_name, tf_n_pos, tf_n_tested, tf_prev) %>%
  mutate(
    country = str_extract(location_name, ",\\s*[^\\(]+") %>%
      str_replace_all("^,\\s*", "")
  )
unique(trachoma_tf$country)

tf_by_age_country <- trachoma_tf %>%
  group_by(year, country, age_years) %>%
  summarise(
    mean_tf_prev = mean(tf_prev, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(tf_by_age_country, aes(x = age_years, y = mean_tf_prev, color = as.factor(year))) +
  geom_line() +
  geom_point() +
  facet_wrap(~ country) +
  labs(
    title = "TF Prevalence by Age and Year",
    x = "Age (years)",
    y = "Mean TF Prevalence (%)",
    color = "Year"
  ) +
  theme_minimal()


## TI POSITIVE
trachoma_ti <- age_bins %>% 
  select(year, age_years, location_name, ti_n_pos, ti_n_tested, ti_prev) %>%
  mutate(
    country = str_extract(location_name, ",\\s*[^\\(]+") %>%
      str_replace_all("^,\\s*", "")
  )
unique(trachoma_ti$country)

ti_by_age_country <- trachoma_ti %>%
  group_by(year, country, age_years) %>%
  summarise(
    mean_ti_prev = mean(ti_prev, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(ti_by_age_country, aes(x = age_years, y = mean_ti_prev, color = as.factor(year))) +
  geom_line() +
  geom_point() +
  facet_wrap(~ country) +
  labs(
    title = "TI Prevalence by Age and Year",
    x = "Age (years)",
    y = "Mean TI Prevalence (%)",
    color = "Year"
  ) +
  theme_minimal()


## PGP3 POSITIVE
trachoma_pgp3 <- age_bins %>% 
  select(year, age_years, location_name, pgp3_n_pos, pgp3_n_tested, pgp3_prev) %>%
  mutate(
    country = str_extract(location_name, ",\\s*[^\\(]+") %>%
      str_replace_all("^,\\s*", "")
  )
unique(trachoma_pgp3$country)

pgp3_by_age_country <- trachoma_pgp3 %>%
  group_by(year, country, age_years) %>%
  summarise(
    mean_pgp3_prev = mean(pgp3_prev, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(pgp3_by_age_country, aes(x = age_years, y = mean_pgp3_prev, color = as.factor(year))) +
  geom_line() +
  geom_point() +
  facet_wrap(~ country) +
  labs(
    title = "PGP3 Prevalence by Age and Year",
    x = "Age (years)",
    y = "Mean PGP3 Prevalence (%)",
    color = "Year"
  ) +
  theme_minimal()


################################################################################
## PROBABILITY ##
################################################################################

data_prob <- age_bins %>%
  select(year, age_years, location_name, pcr_prev, tf_prev, ti_prev, pgp3_prev) %>%
  mutate(
    country = str_extract(location_name, ",\\s*[^\\(]+") %>%
      str_replace_all("^,\\s*", "")
  )
unique(data_prob$country)


prob_age_country <- data_prob %>%
  group_by(year, country, age_years) %>%
  summarise(
    mean_pcr_prev = mean(pcr_prev, na.rm = TRUE),
    mean_tf_prev = mean(tf_prev, na.rm = TRUE),
    mean_ti_prev = mean(ti_prev, na.rm = TRUE),
    mean_pgp3_prev = mean(pgp3_prev, na.rm = TRUE),
    .groups = "drop"
  )


## A focus on Ethiopia for example:
prob_ethopia <- prob_age_country %>%
  filter(str_detect(country, "Ethiopia"))

eth_long <- prob_ethopia %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "marker",
    values_to = "prevalence"
  ) %>%
  filter(!is.na(prevalence))

ggplot(eth_long, aes(x = age_years, y = prevalence, color = marker, group = marker)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ year) +
  labs(
    title = "Trachoma Marker Prevalence in Ethiopia by Age Group and Year",
    x = "Age Group",
    y = "Prevalence (%)",
    color = "Marker"
  ) +
  theme_minimal()
#issues 2017



