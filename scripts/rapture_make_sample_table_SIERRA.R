# RAPTURE bamlist & bamclst creation for different groups/sites
# 2017-Jun

# This script can be used to generate bamlists and bamlist_clst files for different subsampled groups.
# The following pieces should have already occurred and the subsampled bamlists should already be locally 
# available in a "data_output" folder within your R Project.

# STEPS 1A -- 1C HAPPEN ON TERMINAL/COMMAND LINE
# 1A. SUBSAMPLE ON FARM ----------------------------------------------------

# sbatch -p high -t 24:00:00 02b_run_subsample.sh bam_sort_list 30000

# 1B. MAKE BAMLIST_FLT FILE ------------------------------------------------

# ls *flt_100000* > bamlist_flt_100k

# 1C. USE SFTP TO GRAB FILE ------------------------------------------------

# cd to location where your RProject/data_output folder lives
# sftp -P 2022 USERNAME@farm.cse.ucdavis.edu
# cd to location where the subsampled bamlist lieves (e.g., cd projects/rana_rapture/fastq)
# use GET to pull files from cluster/farm to your local drive and PUT for opposite
# get bamlist_flt_100k

# 2. LOAD LIBRARIES ----------------------------------------------------------

suppressMessages({
  library(tidyverse);
  library(lubridate);
  library(magrittr)
})

options(scipen = 12) # to avoid issues with paste functions, & joining on subsample number

# 3. GET DATA ----------------------------------------------------------------

# set subsample number (so need corresponding subsampled bamlist locally, see Steps 1A-1C)
bamNo<-100
# set site name (will be appended into filename)
site<-"RABO_sierra_v3" # sierra_basin, sierra, sierra_v2

# METADATA
metadat<- read_csv("data_output/rapture06_metadata_revised.csv") %>% arrange(Seq)

# somehow dropped some metadat need to fix ecoregions:
# metadat <- metadat %>% 
#   mutate(EcoRegion = case_when(
#     grepl("^BEAR|^MFA|^NFA|^DEER|^MFY|^NFMFA|^SFA|^SFY", Locality) ~ "Sierra Nevada",
#     TRUE ~ EcoRegion
#   ))
# 
# write_csv(metadat, path = "data_output/rapture06_metadata_revised.csv")

# need to change Seq column to match the new merged one:
metadat$Seq <- gsub(pattern = "SOMM165", replacement = "SOMM163", x = metadat$Seq)

# 3b. SUBSAMPLED BAMLIST --------------------------------------------------

# this is the subsampled bamlist: e.g., "bamlist_flt_50k"
# modify per user's naming convention

#bams <- read_tsv(paste0("data_output/bamlists/bamlist_flt_",bamNo,"k"), col_names = F)

bams <- read_tsv(paste0("data_output/bamlists/bamlist_flt_mrg_",bamNo,"k"), col_names = F)

# remove the *000 component for join, requires fixing scipen for digits
subsamp<-bamNo*1000

# bamlist
bams$X1<-gsub(pattern = paste0(".sortflt.mrg_",subsamp,".bam"), replacement = "", bams$X1)

# 4. FILTER BY: -----------------------------------------------------------

# any number of filtering options can occur...below a few different versions

# 4a. FILTER BY HUC OR REGIONS --------------------------------------------

# By HUC8
#dat<- filter(metadat, HU_8_NAME %in% c("Upper Yuba", "South Fork American"))
#dat<- filter(metadat, HU_8_NAME %in% c("North Fork American", "South Fork American") & SPP_ID=="RABO")

# Yuba
# dat <- filter(metadat, SPP_ID=="RABO", HU_10_NAME %in% c("Lower North Yuba River", "Lower South Yuba River", "Middle Yuba River", "Middle North Yuba River", "Deer Creek", "Upper South Yuba River"))

# EcoREgions

# Sierras
# dat <- filter(metadat, SPP_ID=="RABO" | SPP_pc1=="RABO") %>%
#   filter(EcoRegion=="Sierra Nevada")

# Sierras_v2 (minus Stan and Tuo samples)
# dat <- filter(metadat, SPP_ID=="RABO" | SPP_pc1=="RABO") %>%
# filter(EcoRegion=="Sierra Nevada" & !HU_6_NAME=="San Joaquin")

# Sierras_v3 (minus Stan and Tuo samples, and Deer Ck and SFA samples
dat <- filter(metadat, SPP_ID=="RABO" | SPP_pc1=="RABO") %>%
  filter(EcoRegion=="Sierra Nevada" & !HU_6_NAME=="San Joaquin" & !grepl("^DEER-|^SFA-",Locality))

# Sierras/BASIN
# dat <- filter(metadat, SPP_ID=="RABO" | SPP_pc1=="RABO") %>%
#  filter(EcoRegion=="Sierra Nevada" | EcoRegion=="Sierra/Basin Range")

# be aware outliers: 
  ## BEAR-MISS Canyon RAP-278
  ## SFY-SpringCk RAP 122

# Feather
#dat <- filter(metadat, SPP_ID=="RABO" | SPP_pc1=="RABO") %>%
 # filter(EcoRegion=="Sierra/Basin Range")

# North Coast
# dat <- filter(metadat, SPP_ID=="RABO" & 
#                 EcoRegion=="North Coast"| 
#                 EcoRegion=="Klamath/North Coast" |
#                 EcoRegion=="Cascades" |
#                 EcoRegion=="Southern OR Coastal" |
#                 EcoRegion=="Northern CA Coastal Foothills")

# 4e. FILTER BY REG UNREG SIERRAS -----------------------------------------

# filter to RABO only and sites in American, Bear, Feather, Yuba watersheds

#dat <- filter(metadat, SPP_ID=="RABO", grepl('^NFF|NFA|NFY|MFY|MFA|RUB|SFA|SFY|NFMFA|BEAR',Locality))

# 5. JOIN WITH FLT SUBSAMPLE LIST -----------------------------------------

# check col names for join...should be BAMFILE name with plate/well ID
dfout <- inner_join(dat, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# select cols 
dfout <- dfout %>% select(Seq, SampleID, LabID, SPP_ID, YYYY, lat, lon, Locality)
glimpse(dfout)
table(dfout$SPP_ID)

# group and tally
#dfout %>% group_by(Locality) %>% tally %>% write_csv("data_output/rapture_dfout_tally_sierra_V3.csv")

dfout_tally <- dfout %>% group_by(Locality) %>% tally

# now do same with dat:
dat <- dat %>% select(Seq, SampleID, LabID, SPP_ID, YYYY, lat, lon, Locality)
table(dat$SPP_ID)

dat_tally <- dat %>% group_by(Locality) %>% tally

dat_tally

# join
sites_samples <- full_join(dat_tally, dfout_tally, by="Locality")

# change col names
sites_samples <- sites_samples %>% 
  rename(n_initial=n.x, n_retained=n.y) %>% 
  mutate(Locality = gsub("-", "_", Locality),
         Locality = toupper(Locality))

# write out
write_csv(sites_samples, path = "data_output/samples_tally_sierra_v3_raw.csv")

# combine EUCH from EUCH_US, EUCH_DS, NFNFA, Fix RUB_LC_US, delete MFA_US_R, delete MFF site, drop BEAR_MISSOURI canyon, MFY_grizzlyCk, and MFY Remmington,
# drop NFY_slate_onion and combine NFY_slate and NFY_slate_cgrav
# combine SFY fall ck and SFY into SFY
# drop SFY_humbug, and SFY thimb
# drop SFY Scotchman and SFY mckinnon because too few samples

# MAKE TABLE 1 ------------------------------------------------------------

tab1 <- read_csv(file = "data_output/samplesites_Table1.csv")
sites_samples <- read_csv("data_output/samples_tally_sierra_v3_revised.csv") %>% 
  rename(SiteID=Locality)

# append samples to that
tab1_v2 <- full_join(sites_samples, tab1) %>% 
  rename(SiteName=SiteID,
         SiteID=siteID3)

# rearrange
tab1_v2 <- tab1_v2 %>% select(SiteName,SiteID, Locality:NHD_StO,NHD_DA_km2, n_initial, n_retained)

write_csv(tab1_v2, "data_output/samplesites_Table1_v2.csv")

