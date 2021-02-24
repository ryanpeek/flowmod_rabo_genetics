# get years for samples

library(rio)
library(tidyverse)

main_db <- import("data_output/rapture_master_20190720.csv") %>% 
  select(sample_id, yr_collected)

dat_s1 <- import("data_output/S1_all_samples.csv") %>% janitor::clean_names()

# join on sample_id
df_all <- left_join(dat_s1, main_db, by="sample_id") %>% arrange(sample_id)


write_csv(df_all, "manuscript/ECOSPHERE/updated_tableS2_rapture_samples_w_yr.csv")


# GET A BREAKDOWN OF YEARS ------------------------------------------------


dat <- import("manuscript/ECOSPHERE/updated_tableS2_rapture_samples_w_yr.csv")

prop2014 <- dat %>% filter(yr_collected>2014) %>% count() / count(dat)

