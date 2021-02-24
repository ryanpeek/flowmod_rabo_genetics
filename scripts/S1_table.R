# make supplemental tables:

library(tidyverse)

# read in hbs samples
hbs <- read_csv("data/Samples_HBS_airtable_20161011.csv")

# read in s1
load("data_output/S1_table_old.rda")

# rename
s1 <- supplTAB
s1 <- select(s1, -reads_mean)
rm(supplTAB)



# Join Tables -------------------------------------------------------------

# join any samples from HBS as attrib.
names(hbs)

# join v1
hbs_sel <- select(hbs, SampleID, HBS_ID, Collection, Tissue) %>% 
  mutate(Tissue = "tissue") %>% rename(SampleType = Tissue)
full_s1 <- left_join(s1, hbs_sel) 

table(full_s1$Collection)

# delete RAP1649
full_s1 <- full_s1 %>% filter(!SampleID=="RAP1649")

glimpse(full_s1)

# add RAP to other Collection
full_s1 <- full_s1 %>% 
  mutate(Collection = case_when(
    grepl("^RAP1", SampleID) ~ "CSG", # Caren Goldberg
    grepl("RAP6|RAP7|RAP8|RAP9|AA", SampleID) ~ "RAP",
    TRUE ~ Collection)
  )

# change RAP1649 to "NFY-SLATE-CGRav"
full_s1 <- full_s1 %>% mutate(Locality = replace(Locality, Locality=="NFY-SLATE", "NFY-SLATE-CGRav"))


# tissue in HBS: 
full_s1 <- full_s1 %>% 
  mutate(SampleType = case_when(
    grepl("^RAP1", SampleID) ~ "tissue", # Caren Goldberg
    TRUE ~ SampleType))


# join with metadat for tissue/buccal?
rap_samples <- read_csv("data/Samples_RAP_airtable_20170310.csv") %>% 
  select(SampleID, Swab) %>% 
  mutate(sampletype = case_when(
    grepl("B", Swab) ~ "buccal",
    TRUE ~ "tissue")
  )

# join back:
full_s1_rev <- left_join(full_s1, rap_samples) %>% 
  mutate(SampleType=if_else(is.na(SampleType), sampletype, SampleType)) %>% 
  select(-Swab, -sampletype)


#write_csv(full_s1_rev, "data_output/S1_all_samples.csv")


table(full_s1_rev$SampleType)

table(full_s1_rev$Collection)


full_s1_rev %>% group_by(Locality) %>% tally %>% print(n=Inf) 



