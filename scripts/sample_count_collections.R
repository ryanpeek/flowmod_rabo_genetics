# site and sample assessment:

# how many and who?

library(tidyverse)


# READ IN DATA ------------------------------------------------------------
df <- read_csv("data_output/rapture06_metadata_revised.csv")
head(df)

# total 1444 samples
# RABO 1103 (all)

plates <- read_csv("data/Plates_RADseq-RAPTURE.csv")
head(plates)


# MERGE IN PLATE COLLECTION INFO ------------------------------------------

df2 <- left_join(df, plates[,c(2:5)], by=c("PlateID", "RAD_ID")) %>% select(-Collection.x) %>% 
  rename(Collection=Collection.y)

table(df2$Collection)

# COUNT & TALLY  ----------------------------------------------------------

df_rb <- filter(df2, SPP_ID=="RABO", !grepl("WSU", Collection))

# RABO 985  (RABO_all_50k used 809 individuals, so ~80+ percent sample success?)

## Sierras Specifically:
table(df_rb$EcoRegion)
df_rb %>% filter(EcoRegion=="Sierra Nevada") %>% tally # 312
