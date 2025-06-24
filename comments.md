
Can you document where you got these data?
trachoma_sero_transmission_analysis_cluster.csv
trachoma_sero_transmission_analysis_indiv.csv
trachoma_sero_transmission_analysis_study.csv

Please avoid long paths (and paths with spaces) in your scripts. 
* cluster_data <- read_csv("Documents/AIMS/MMED/MMED R/MMED TRACHOMA/…

It's generally good to:
* Control where your script is running
* Put the data somewhere nearby so you can access with local paths 
e.g., `cluster_data <- read_csv("data/trachoma_sero_transmission_analysis_cluster.csv")`
* Tell git to ignore the data path

I'm also a big fan of breaking scripts in pieces. Once you have a script that reads and cleans the data, save that data, and start a new script that reads the data. It's a good idea to often run scripts from beginning to end.

It's good to comment only when it will help the flow, or when something is deep: 

```
summary(cluster_data)   # Summarize cluster-level data
mutate(age = as.numeric(age_years)) %>%        # Ensure age is numeric
```

are examples where the comments aren't necessary.

This script seems to get very deep into a particular approach, maybe before questioning enough whether it's the best approach to take? Maybe? I'm not 100% following the boot stuff at the end.
