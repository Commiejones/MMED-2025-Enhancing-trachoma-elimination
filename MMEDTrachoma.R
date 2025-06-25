library(readr)       
library(dplyr)       
library(tidyr)       
library(purrr)       
library(ggplot2) 

cluster_data <- read_csv("trachoma_sero_transmission_analysis_cluster.csv", show_col_types = F)
individual_data <- read_csv("trachoma_sero_transmission_analysis_indiv.csv", show_col_types = F)
study_data <- read_csv("trachoma_sero_transmission_analysis_study.csv", show_col_types = F)

##join the cluster data by cluster_id to individual data and keep the following columns
#cluster_id,pcr_n_tested,pcr_n_pos,tf_n_tested,tf_n_pos,,pgp3_n_tested,pgp3_n_pos, age_months, age_years,
#location_name, year, 
cluster_data1 <- cluster_data %>%
  select(cluster_id, pcr_n_tested, pcr_n_pos, tf_n_tested, tf_n_pos, pgp3_n_tested, 
         pgp3_n_pos)
individual_data1 <- individual_data %>%
  select(cluster_id, age_years, location_name, year) %>%
  mutate(location_name = as.character(location_name))

#join the two data frames by cluster_id
joined_data <- individual_data1 %>%
  left_join(cluster_data1, by = "cluster_id") %>%
  mutate(location_name = as.character(location_name))

#sorted by age and year
joined_data <- joined_data %>%
  arrange(age_years, year)

#put age structure 1-5 and 6-9 together

joined_data <- joined_data %>%
  mutate(age_years = ifelse(age_years <= 5, "1-5", "6-9")) %>%
  group_by(location_name, year, age_years) %>%
  summarise(
    across(starts_with("pcr_n_"), sum, na.rm = TRUE),
    across(starts_with("tf_n_"), sum, na.rm = TRUE),
    across(starts_with("pgp3_n_"), sum, na.rm = TRUE),
    .groups = 'drop'
  )



#Extract the country from the location_name 
joined_data <- joined_data %>%
  mutate(country = sub(".*,\\s*", "", location_name)) %>%
  mutate(country = sub("\\s*\\(.*?\\)", "", country)) %>%
  mutate(country = as.character(country))


#rename columns 
joined_data <- joined_data %>%
  rename(
    `PCR total` = pcr_n_tested,
    `PCR pos` = pcr_n_pos,
    `TF total` = tf_n_tested,
    `TF pos` = tf_n_pos,
    `PGP3(Serology) total` = pgp3_n_tested,
    `Sero pos` = pgp3_n_pos
  )

#calculate prevalence

df <- joined_data %>%
  mutate(
    TF_prev = `TF pos` / `TF total`,
    Sero_prev = `Sero pos` / `PGP3(Serology) total`,
    PCR_prev = `PCR pos` / `PCR total`
  )
summary_by_location_age <- df %>%
  group_by(year,location_name, age_years) %>%
  summarise(
    mean_TF_prev = mean(TF_prev, na.rm = TRUE),
    mean_Sero_prev = mean(Sero_prev, na.rm = TRUE),
    mean_PCR_prev = mean(PCR_prev, na.rm = TRUE)
  )

summary_by_country_age <- df %>%
  group_by(country, location_name, age_years,year) %>%
  summarise(
    mean_TF_prev = mean(TF_prev, na.rm = TRUE),
    mean_Sero_prev = mean(Sero_prev, na.rm = TRUE),
    mean_PCR_prev = mean(PCR_prev, na.rm = TRUE)
  )


# scatter plot that shows the trend for prevalence on the y axis and color by TF, PCR and Sero prevalence for each age group
ggplot(df, aes(x = year)) +
  geom_point(aes(y = TF_prev, color = "TF Prevalence")) +
  geom_point(aes(y = Sero_prev, color = "Seroprevalence")) +
  geom_point(aes(y = PCR_prev, color = "PCR Prevalence")) +
  labs(title = "Prevalence Trends by Year",
       x = "Age Group", y = "Prevalence") +
  scale_color_manual(values = c("TF Prevalence" = "blue", 
                                "Seroprevalence" = "red", 
                                "PCR Prevalence" = "green")) +
  theme_minimal()

# scatter plot that shows the trend for prevalence on the y axis and color by TF, PCR and Sero prevalence for each age group and 
#country and year

all_years <- sort(unique(df$year))
all_years
for (y in all_years) {
  plt <- ggplot(df %>% filter(year == y), aes(x = age_years)) +
    geom_point(aes(y = TF_prev, color = "TF Prevalence")) +
    geom_point(aes(y = Sero_prev, color = "Seroprevalence")) +
    geom_point(aes(y = PCR_prev, color = "PCR Prevalence")) +
    facet_wrap(~ location_name) +
    labs(title = paste("Prevalence Trends by Age Group and Country in", y),
         x = "Age Group", y = "Prevalence") +
    scale_color_manual(values = c("TF Prevalence" = "blue", 
                                  "Seroprevalence" = "red", 
                                  "PCR Prevalence" = "green")) +
    theme_minimal()
  print(plt)
}

