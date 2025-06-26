library(shellpipes)
library(dplyr)

dat <- rdsRead()
summary(dat)

loc <- (dat
	|> group_by(location_name)
	|> summarise(
		pcrN = sum(!is.na(pcr))
		, pcr = mean(pcr, na.rm=TRUE)
		, tfN = sum(!is.na(tf))
		, tf = mean(tf, na.rm=TRUE) ## , tf, ti, abx = pgp3_pos
		, abxN = sum(!is.na(abx))
		, abx = mean(abx, na.rm=TRUE) ## , tf, ti, abx = pgp3_pos
	)
	|> ungroup()
)

locYear <- (dat
	|> group_by(location_name, year)
	|> summarise(
		pcrN = sum(!is.na(pcr))
		, pcr = mean(pcr, na.rm=TRUE)
		, tfN = sum(!is.na(tf))
		, tf = mean(tf, na.rm=TRUE) ## , tf, ti, abx = pgp3_pos
		, abxN = sum(!is.na(abx))
		, abx = mean(abx, na.rm=TRUE) ## , tf, ti, abx = pgp3_pos
	)
	|> ungroup()
)

rdaSave(loc, locYear)

