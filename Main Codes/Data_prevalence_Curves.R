################################################################################
## PACKAGES AND LIBRARIES ##
################################################################################
library(readr)       
library(dplyr)       
library(tidyr)       
library(purrr)       
library(ggplot2)
library(stringr)
library(scales)

################################################################################
## LOAD DATA ##
################################################################################
indiv_data <- read_csv("trachoma_sero_transmission_analysis_indiv.csv")

################################################################################
## MODEL - Wag Hemra, Ethiopia (WUHA) ##
################################################################################

years <- 2013:2019

for (yr in years) {
  
  # Summarise data by age
  prev_by_age <- indiv_data %>%
    filter(
      year == yr,
      location_name == 'Wag Hemra, Ethiopia (WUHA)',
      age_years %in% 1:9
    ) %>%
    drop_na(tf, pgp3_pos) %>%
    group_by(age_years) %>%
    summarise(
      total = n(),
      positive_sero = sum(pgp3_pos == 1),
      positive_tf = sum(tf == 1),
      prevalence_sero = positive_sero / total,
      prevalence_tf = positive_tf / total
    ) %>%
    arrange(age_years)
  
  # Skip year if there's no data for the year
  if (nrow(prev_by_age) == 0) next
  
  #Long format
  prev_long <- prev_by_age %>%
    select(age_years, prevalence_sero, prevalence_tf) %>%
    pivot_longer(
      cols = starts_with("prevalence"),
      names_to = "indicator",
      values_to = "prevalence"
    ) %>%
    mutate(
      indicator = recode(indicator,
                         prevalence_sero = "Seroprevalence (Pgp3)",
                         prevalence_tf = "TF Prevalence")
    )
  
  #Plot
  p <- ggplot(prev_long, aes(x = age_years, y = prevalence, color = indicator)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = 1:9) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_color_manual(
      values = c(
        "Seroprevalence (Pgp3)" = "black",
        "TF Prevalence" = "red"
      )
    ) +
    labs(
      title = paste0("Age-Specific Seroprevalence and TF Prevalence in WUHA (", yr, ")"),
      subtitle = "Children Aged 1–9 in Wag Hemra, Ethiopia",
      x = "Age (Years)",
      y = "Prevalence (%)",
      color = "Indicator"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  print(p) 
  
  ggsave(
    filename = paste0("WUHA_prevalence_", yr, ".jpeg"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}



################################################################################
## MODEL - Kongwa, Tanzania ##
################################################################################

years <- 2013:2019

for (yr in years) {
  
  #Summarise data by age
  prev_by_age <- indiv_data %>%
    filter(
      year == yr,
      location_name == 'Kongwa, Tanzania',
      age_years %in% 1:9
    ) %>%
    drop_na(tf, pgp3_pos) %>%
    group_by(age_years) %>%
    summarise(
      total = n(),
      positive_sero = sum(pgp3_pos == 1),
      positive_tf = sum(tf == 1),
      prevalence_sero = positive_sero / total,
      prevalence_tf = positive_tf / total
    ) %>%
    arrange(age_years)
  
  # Skip year if there's no data for the year
  if (nrow(prev_by_age) == 0) next
  
  # Long format
  prev_long <- prev_by_age %>%
    select(age_years, prevalence_sero, prevalence_tf) %>%
    pivot_longer(
      cols = starts_with("prevalence"),
      names_to = "indicator",
      values_to = "prevalence"
    ) %>%
    mutate(
      indicator = recode(indicator,
                         prevalence_sero = "Seroprevalence (Pgp3)",
                         prevalence_tf = "TF Prevalence")
    )
  
  #Plot
  p <- ggplot(prev_long, aes(x = age_years, y = prevalence, color = indicator)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = 1:9) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_color_manual(
      values = c(
        "Seroprevalence (Pgp3)" = "black",
        "TF Prevalence" = "red"
      )
    ) +
    labs(
      title = paste0("Age-Specific Seroprevalence and TF Prevalence in Kongwa, Tanzania (", yr, ")"),
      subtitle = "Children Aged 1–9 in Kongwa, Tanzania",
      x = "Age (Years)",
      y = "Prevalence (%)",
      color = "Indicator"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)  
  
  ggsave(
    filename = paste0("Tanzania_prevalence_", yr, ".jpeg"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
}
