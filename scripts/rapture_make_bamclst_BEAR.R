# RAPTURE bamlist & bamclst creation for different POPS/SUBPOPS
# This script can be used to generate bamlists and bamlist_clst files for different subsampled groups.
# The following pieces should have already occurred and the subsampled bamlists should already be locally 
# available in a "data_output" folder within your R Project.

# Updated: Thu Sep 21 20:36:38 2017 ------------------------------

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
bamNo<-50

# set site names (will be appended into filename)
site<-"RABO_bear"

# METADATA
metadat<- read_csv("data_output/rapture06_metadata_revised.csv") %>% arrange(Seq)

# need to change Seq column to match the new merged one:
metadat$Seq <- gsub(pattern = "SOMM165", replacement = "SOMM163", x = metadat$Seq)

# 3b. SUBSAMPLED BAMLIST --------------------------------------------------

bams <- read_tsv(paste0("data_output/bamlists/bamlist_flt_mrg_",bamNo,"k"), col_names = F)
subsamp<-bamNo*1000 # remove the *000 component for join
bams$X1<-gsub(pattern = paste0(".sortflt.mrg_",subsamp,".bam"), replacement = "", bams$X1)

# 4. FILTER BY SITE NAME -------------------------------------------------

# POP: this is for main bamlist
dat <- filter(metadat, SPP_ID=="RABO", grepl('^BEAR', Locality), SampleID!="RAP-278") 

# SUBpop: Splits by Locality col into list of unique dataframes
dat2 <- split(x = dat, dat$Locality)
names(dat2)<-gsub(pattern = "-", x = tolower(names(dat2)), replacement = "_")

# split apart into indiv data frames in GlobalEnv
list2env(x = dat2, .GlobalEnv)

# Split/Rename/Combine here if necessary
bear_stho <- bear_sth2
rm(bear_sth2)
bear_stho <- bind_rows(bear_stho, bear_sthc)
rm(bear_sthc)
bear_stho_haw <- bear_stha
rm(bear_stha)
bear_sthohaw <- bind_rows(bear_stho, bear_stho_haw)
rm(bear_miss_cnyn)

# see how many groups:
#dat %>% group_by(Locality) %>% tally() %>% as.data.frame()


# rm all the data.frames starting with site
#rm(list = ls(pattern = "mfa"))

# 5a. JOIN & MAKE BAMLIST/CLST FILES: POPS --------------------------------

# WRITE BAMLIST
dfout <- inner_join(dat, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# Write to bamlist for angsd call:
write_delim(as.data.frame(paste0(dfout$Seq, ".sortflt.mrg_",subsamp,".bam")),
            path = paste0("data_output/bamlists/bamlist_mrg_",site,"_",bamNo,"k"), col_names = F)

names(dfout) # check and see what you've got

# WRITE BAMCLST

# FINE SCALE
clst_out<-dfout %>%
  dplyr::select(HU_10_NAME, SampleID, Locality) %>%
  dplyr::rename(FID=HU_10_NAME, IID=SampleID, CLUSTER=Locality)
length(unique(clst_out$CLUSTER))
length(unique(clst_out$FID))

# MED SCALE
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

write_delim(clst_out, path=paste0("data_output/bamlists/bamlist_mrg_",site,"_",bamNo,"k_clst"))

# 5b. JOIN & MAKE BAMLIST/CLST FILES: SUBPOPS -----------------------------

# FOR SUBPOPS
sites<-ls(pattern = "bear") # need to make a list of site names

# function for subpops:
make_bamclst <- function(sitename) {
  dfout <- inner_join(get(sitename), bams, by=c("Seq"="X1")) %>% arrange(Seq)
  
  # Write to bamlist for angsd call:
  write_delim(as.data.frame(paste0(dfout$Seq, ".sortflt.mrg_",subsamp,".bam")),
              path = paste0("data_output/bamlists/bamlist_mrg_",
                            sitename,"_",bamNo,"k"), col_names = F)
  clst_out<-dfout %>%
    dplyr::select(HUC_10, SampleID, Locality) %>%
    dplyr::rename(FID=HUC_10, IID=SampleID, CLUSTER=Locality)
  head(clst_out)
  write_delim(clst_out, 
              path=paste0("data_output/bamlists/bamlist_mrg_",
                          sitename,"_",bamNo,"k_clst"))
}

map(sites, make_bamclst)

# 6. TERMINAL SFTP --------------------------------------------------------

# cd Documents/github/rabo_regulation/data_output/bamlists
# farmer #(sftp)
# cd projects/rana_rapture/MERGED/bamlists/subpops
sitename <- "RABO_bear"
paste0("put bamlist_mrg_",sitename,"*",bamNo,"k*")

# cd projects/rana_rapture/MERGED/bamlists/POPS
#paste0("put bamlist_mrg_",site,"_*",bamNo,"k*")

# 7. BASH: PCA CALC SITES ------------------------------------------------

# Use angsd to run pca_calc_sites script: 
lsite<- tolower(site)
paste0("sbatch -p high -t 24:00:00 03_pca_calc_sites.sh bamlist_mrg_",site,"_",bamNo,"k", " ", lsite, "_",bamNo,"k", " bamlists/POPS")

# 8a. BASH: Thetas/SFS ---------------------------------------------------

## in bamlists/subpops:
# ls bamlist_mrg_bear* | grep "k$" | sed 's/bamlist_mrg_//g' > thetalist_bear

## mv list out:
# mv thetalist_bear ../../

# sbatch -p high -t 12:00:00 06_theta_sfs.sh thetalist_bear

# 8. BASH: Calc Pairwise FST ---------------------------------------------

sites<-ls(pattern = "bear") 

# get all unique pairs:
sitepairs <- combn(x = sites, m = 2)

for(i in 1:ncol(sitepairs)){
  cat(paste0("sbatch -p high -t 12:00:00 07_get_fst.sh ",sitepairs[,i][1],"_",bamNo,"k", " ",sitepairs[,i][2],"_",bamNo,"k"))
  cat("\n")
}

# can check datestamp & sort by time: ls -lt FILE*

sbatch -p high -t 12:00:00 07_get_fst.sh bear_100k bear_grhn_100k
sbatch -p high -t 12:00:00 07_get_fst.sh bear_100k bear_stho_100k
sbatch -p high -t 12:00:00 07_get_fst.sh bear_100k bear_stho_haw_100k
sbatch -p high -t 12:00:00 07_get_fst.sh bear_100k bear_sthohaw_100k
sbatch -p high -t 12:00:00 07_get_fst.sh bear_grhn_100k bear_stho_100k
sbatch -p high -t 12:00:00 07_get_fst.sh bear_grhn_100k bear_stho_haw_100k
sbatch -p high -t 12:00:00 07_get_fst.sh bear_grhn_100k bear_sthohaw_100k
sbatch -p high -t 12:00:00 07_get_fst.sh bear_stho_100k bear_stho_haw_100k


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




