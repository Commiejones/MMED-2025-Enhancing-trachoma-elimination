library(readr)       
library(dplyr)       

library(shellpipes)

dat <- csvRead()

dat |> mutate_if(is.character, as.factor) |> summary() |> print()

dat <- (dat
	|> transmute(location_name, year
		, household_id, individual_id
		, pcr, tf, ti 
		, age = ifelse(is.na(age_months), 12*age_years+6, age_months)
	)
)
