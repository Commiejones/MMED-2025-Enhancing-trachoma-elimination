library(dplyr)
library(readr)

dat <- read_csv("data/trachoma_sero_transmission_analysis_indiv.csv")

print(dat
	|> mutate_if(is.character, as.factor)
	|> summary()
)

## There are 22 combinations of sites and years
print(dat
	|> select(study_id, location_name, location_year_name, year)
	|> distinct()
)

## These are captured by location_name and year
print(dat
	|> select(location_name, year)
	|> distinct()
)

## There are 14 overall combinations of sites
## There is a glitch here: Kongwa corresponds to two different studies
## while TCC-Ethiopia2017 corresponds to four different locations
## I suggest indexing by location_name and year only
## We can ignore for now that Kongwa is two different studies
print(dat
	|> select(study_id, location_name)
	|> distinct()
)
