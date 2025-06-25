library(readr)       
library(dplyr)       

library(shellpipes)

dat <- csvRead()

dat |> mutate_if(is.character, as.factor) |> summary() |> print()

dat <- (dat
	|> select(location_name, year
		, household_id, individual_id
		, age_months, age_years
		, pcr, tf, ti, pgp3_mfi, pgp3_mfi_nonneg, pgp3_mfi_log10, pgp3_elisa, pgp3_mfi_cutoff, pgp3_pos, pgp3_minobs 

)
