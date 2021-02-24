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


# FOR SAMTOOLS EVAL -------------------------------------------------------

samsite <- "samtools_sierra_v3"
#dat2 <- dat %>% group_by(Locality) %>% filter(n() > 2) # this doesn't work because of site combos...need to use specific filter:

dat %>% group_by(Locality) %>% tally() %>% print(n=Inf)

# drop: BEAR-MISS_CNYN, MFA-US-R, MFF-SMTCk, MFY-GRIZZLY_CK, MFY-Remmington, NFY-SLATE-Onion, SFY-FallCk, SFY-HUMBUG, SFY-MCKI, SFY-THIMC, SFY-Scotchman

dat2 <- dat %>% filter(!Locality %in% c("BEAR-MISS_CNYN", "MFA-US-R", "MFF-SMTCk", "MFY-GRIZZLY_CK", "MFY-Remmington", "NFY-SLATE-Onion", "SFY-FallCk", "SFY-HUMBUG", "SFY-MCKI", "SFY-THIMC", "SFY-Scotchman"))

dat2 %>% group_by(Locality) %>% tally() %>% print(n=Inf)

# Write to bamlist for samtools (orig unfiltered):
write_delim(as.data.frame(paste0(dat2$Seq, ".sortflt.mrg_",subsamp,".bam")),
            path = paste0("data_output/bamlists/bamlist_mrg_",samsite,"_",bamNo,"k"), col_names = F)

write_delim(as.data.frame(paste0(dat2$Seq, ".sortflt.mrg.bam")),
            path = paste0("data_output/bamlists/bamlist_mrg_",samsite,"_",bamNo,"k_thresh"), col_names = F)


# 5. JOIN WITH FLT SUBSAMPLE LIST -----------------------------------------

# check col names for join...should be BAMFILE name with plate/well ID
dfout <- inner_join(dat, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# 6. WRITE TO BAMLIST -----------------------------------------------------

# Write to bamlist for angsd call:
write_delim(as.data.frame(paste0(dfout$Seq, ".sortflt.mrg_",subsamp,".bam")),
            path = paste0("data_output/bamlists/bamlist_mrg_",site,"_",bamNo,"k"), col_names = F)

# Write to bamlist_thresh for angsd PCA IBS call (remove subsample_append:
write_delim(as.data.frame(paste0(dfout$Seq, ".sortflt.mrg.bam")),
            path = paste0("data_output/bamlists/bamlist_mrg_",site,"_",bamNo,"k_thresh"), col_names = F)


# run the run cnts file:
"srun -t 12:00:00 01d_run_counts.sh bamlists/POPS/bamlist_mrg_samtools_sierra_v3_100k_thresh
"
# then merge
"paste bamlist_mrg_samtools_sierra_v3_100k_thresh count1_reads.txt count2_mapped.txt count3_paired.txt count4_unique.txt > sierra_samples_v3_fq_counts.txt"

# 7. WRITE CLST FILE ---------------------------------------------------------

names(dfout) # check and see what you've got

# write to file for later use/filtering
#write_rds(x = dfout, path = "data_output/bamlist_dat_RABO_reg_unreg.rds")

# FINE SCALE
clst_out<-dfout %>%
  dplyr::select(HU_10_NAME, SampleID, Locality) %>%
  dplyr::rename(FID=HU_10_NAME, IID=SampleID, CLUSTER=Locality)
length(unique(clst_out$CLUSTER))
length(unique(clst_out$FID))

# clst_out<-dfout %>%
#   dplyr::select(Locality, SampleID, HU_8_NAME) %>%
#   dplyr::rename(FID=Locality, IID=SampleID, CLUSTER=HU_8_NAME)
# length(unique(clst_out$CLUSTER))
# length(unique(clst_out$FID))

# COARSE SCALE
# clst_out<-dfout %>%
#   dplyr::select(HU_8_NAME, SampleID, HUC_6) %>%
#   dplyr::rename(FID=HU_8_NAME, IID=SampleID, CLUSTER=HUC_6)
# length(unique(clst_out$CLUSTER))
# length(unique(clst_out$FID))

# count NAs in your clst file
clst_out %>% filter(is.na(FID)) %>% tally
clst_out %>% filter(is.na(CLUSTER)) %>% tally

# fill NAs with FEATHER
# clst_out <- clst_out %>% 
#   mutate(FID = case_when(
#     is.na(FID) ~ "EBNF Feather",
#     TRUE ~ FID)
#   )

write_delim(clst_out, path=paste0("data_output/bamlists/bamlist_mrg_",site,"_",bamNo,"k_clst"))

# 8. TERMINAL SFTP --------------------------------------------------------

# cd Documents/github/rabo_regulation/data_output/bamlists/
# farmer (sftp)
# cd projects/rana_rapture/MERGED/bamlists/POPS
paste0("put bamlist_mrg_",site,"_*",bamNo,"k*")

# 9. BASH: PCA CALC SITES --------------------------------------------------

# Use angsd to run pca_calc_sites script: 

# create command:
lsite<- tolower(site)

paste0("sbatch -p high -t 24:00:00 03_pca_calc_sites.sh bamlist_mrg_",site,"_",bamNo,"k", " ", lsite, "_",bamNo,"k", " bamlists/POPS")


# 9a. BASH: PCA NEW -------------------------------------------------------

# make new bamlist:
paste0("less bamlist_mrg_",site,"_",bamNo,"k ") # part 1

# part 2 (won't cat or past0)
"| sed 's/\_100000//g' > " 

paste0("bamlist_mrg_",site,"_",bamNo,"k_thresh") # part 3

# need to make new bamlist of all original bams that meet 100k threshold, but aren't subsampled:
# less bamlist_mrg_RABO_sierra_100k | sed 's/\_100000//g' > bamlist_mrg_RABO_sierra_100k_thresh

# check and see if lots of folks running on bigmemm, if so use high.

paste0("sbatch -p high --mem=64G -t 48:00:00 --mail-type ALL --mail-user rapeek@ucdavis.edu 03a_pca_new.sh bamlist_mrg_",site,"_",bamNo,"k_thresh", " ", lsite, "_",bamNo,"k_thresh", " bamlists/POPS")


# 10. BASH: MAKE PCA PLOTS --------------------------------------------------

# cd projects/rana_rapture/fastq
## (sbatch or srun)

lsite<- tolower(site)
paste0("srun -p med -t 24:00:00 04_pca_plot.sh ", lsite, "_",bamNo,"k", " ","bamlist_",site,"_",bamNo,"k_clst")

# sbatch -p med -t 12:00:00 04_pca_plot.sh RABO_nfa_75k bamlist_RABO_nfa_75k_clst

# can check datestamp & sort by time: ls -lt FILE*


# 11. SFTP and GET PDFS ---------------------------------------------------

# then sftp and "get" to get pdfs back to your computer (see step 1C).

# cd Documents/github/rana_rapture/data_output/bamlists
# farmer (sftp)
# cd projects/rana_rapture/fastq/pca_pdfs

# get rabo_nfa_75k* # (this goes from cluster to local)
paste0(lsite,"_",bamNo,"k*")


# MAPS: STATIC ----------------------------------------------------

# map_data("county") %>% 
#   filter(region %in% c("california","oregon")) %>%
#   ggplot() +
#   geom_polygon(aes(x = long, y = lat, group = group)) + 
#   coord_map("gilbert") +
#   geom_point(data=dfout, aes(x=lon, y=lat, fill=HU_6_NAME), size=4, pch=21)+
#   scale_fill_viridis_d("HUC") + 
#   ggtitle(paste0(site, ": ", bamNo,"k (n=",nrow(dfout),")"))

# MAPS: LEAFLET ----------------------------------------

## if you have X/Y for data, look at a map of samples from subsample

# library(sf)
# library(leaflet)
# 
# # hydrology (HUCS)
# h8 <- st_read("data_output/shps/HUC8_named_westcoast.shp", quiet = T)
# h8 <- st_transform(h8, crs = 4326)
# h6 <- st_read("data_output/shps/WBD_HUC6.shp", quiet = T)
# h6 <- st_transform(h6, crs = 4326)
# 
# # make map
# subMAP <- leaflet() %>%
#   addTiles(group = "Basemap") %>%
#   #setView(lng = -120.8, lat = 39, zoom = 5) %>%  # set to Auburn/Colfax, zoom 5 for CA
#   addProviderTiles("Stamen.TopOSMFeatures", group = "OSM Features") %>%
#   addProviderTiles("Esri.WorldTopoMap", group = "Topo") %>%
#   addProviderTiles("Esri.WorldImagery", group = "ESRI Aerial") %>%
# 
#   # huc6
#   addPolygons(data=h6, group="HUC6", color="#473C8B", weight = 1.5,
#               fillColor = "transparent", label = ~paste0(HU_6_Name, ": " ,HUC_6)) %>%
#   hideGroup("HUC6") %>%
#   # huc8
#   addPolygons(data=h8, group="HUC8", color="darkblue", weight = 1.3,
#               fillColor = "transparent", label = ~HU_8_NAME) %>%
#   hideGroup("HUC8") %>%
# 
#   # add scale bar
#   addMeasure(position = "topright",
#              primaryLengthUnit = "kilometers",
#              primaryAreaUnit = "hectares",
#              secondaryAreaUnit = "sqmiles",
#              activeColor = "#3D535D",
#              completedColor = "#7D4479") %>%
# 
#   # add fancy zoom button to recenter
#   # addEasyButton(easyButton(
#   #   icon="fa-globe", title="Zoom to Level 5",
#   #   onClick=JS("function(btn, map){ map.setZoom(5); }"))) %>%
# 
#   # add subsamples
#   addCircleMarkers(data=dfout, group="Samples",
#                    clusterOptions = markerClusterOptions(),
#                    lng = ~lon, lat=~lat,opacity = 0.5,
#                    popup=paste0("Locality: ", dfout$Locality, "<br>",
#                                 "HUC12Name: ", dfout$HU_12_NAME,
#                                 "<br>","SampleID: ",dfout$SampleID,
#                                 "<br>", "SPP_ID: ",dfout$SPP_ID,
#                                 "<br>", "Elev (m): ", dfout$elev_m, "<br>",
#                                 "StreamName: ", dfout$GNIS_NAME),
#                    weight=0.6,radius=10, stroke=TRUE,
#                    fillColor = ifelse(dfout$SPP_ID=="RABO" |
#                                         dfout$SPP_pc1=="RABO", "yellow", "dodgerblue")) %>%
# 
#   # add layer/legend control
#   addLayersControl(
#     baseGroups = c("Basemap", "Topo", "ESRI Aerial", "OSM Features"),
#     overlayGroups = c("HUC6", "HUC8", "Samples"),
#     options = layersControlOptions(collapsed = T))
# 
# subMAP




