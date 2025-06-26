library(shellpipes)
library(ggplot2); theme_set(theme_minimal(base_size=14))
pdf(width=9)

loadEnvironments()

base <- (ggplot()
	+ aes(color=location_name)
	+ geom_point()
	+ geom_line()
	+ geom_abline(slope=1)
)

print(base %+% loc + aes(tf, abx))
print(base %+% loc + aes(pcr, abx))
print(base %+% loc + aes(pcr, tf))

print(base %+% locYear + aes(tf, abx))
print(base %+% locYear + aes(pcr, abx))
print(base %+% locYear + aes(pcr, tf))
