## Conclusion: pgp3_pos is patched together from elisa, pgp3_mfi and the cutoff; not replicating for now (even though I would like to do a ratio thing).
library(readr)       
library(dplyr)       

library(shellpipes)

dat <- csvRead()

dat <- (csvRead()
	|> select(
		pgp3_mfi, pgp3_mfi_nonneg, pgp3_mfi_log10
		, pgp3_elisa, pgp3_mfi_cutoff, pgp3_pos, pgp3_minobs 
	)
)

summary(dat)

quit()

(dat |> transmute(
	pgp3_pos
	, comp = as.numeric(pgp3_mfi>pgp3_mfi_cutoff)
	, diff = comp-pgp3_pos
)) |> summary()
