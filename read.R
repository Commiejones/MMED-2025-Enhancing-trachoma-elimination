library(readr)       
library(dplyr)       

library(shellpipes)

dat <- csvRead()

dat |> mutate_if(is.character, as.factor) |> summary() |> print()

dat <- (dat
)
study_id, cluster_id, household_id, individual_id, survey, year, location_name, location_year_name, mda, age_months, age_years, pcr, tf, ti, pgp3_mfi, pgp3_mfi_nonneg, pgp3_mfi_log10, pgp3_elisa, pgp3_mfi_cutoff, pgp3_pos, pgp3_minobs 
