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
bamNo<-100

# set site names (will be appended into filename)
site<-"RABO_yuba_subsamp"

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
dat <- filter(metadat, SPP_ID=="RABO", grepl('^SFY|^MFY|^NFY|^DEER', Locality), SampleID!="RAP-122") 
#dat <- filter(metadat, SPP_ID=="RABO", grepl('^NFY', Locality), SampleID!="RAP-122") 

# see how many groups:
#dat %>% group_by(Locality) %>% tally() %>% as.data.frame()

# SUBpop: Splits by Locality col into list of unique dataframes
dat2 <- split(x = dat, dat$Locality)
names(dat2)<-gsub(pattern = "-", x = tolower(names(dat2)), replacement = "_")

# split apart into indiv data frames in GlobalEnv
list2env(x = dat2, .GlobalEnv)

# Split/Rename/Combine here if necessary

# need to drop outlier from SFY-SPRING (RAP-122)
sfy_springck <- filter(sfy_springck, SampleID!="RAP-122")

# create filtered version based on PCA but retain full group as well
mfy_oregck_flt <- filter(mfy_oregck, SampleID!="RAP-149", 
                         SampleID!="RAP-184", SampleID!="RAP-201",
                         SampleID!="RAP-244", SampleID!="RAP-247", 
                         SampleID!="RAP-248")

# combine slate samples
nfy_slate <- bind_rows(nfy_slate, nfy_slate_cgrav, nfy_slate_onion)
rm(nfy_slate_cgrav,nfy_slate_onion)

# combine thimble and logan, all part of LOGAN canyon
sfy_loga <- bind_rows(sfy_thimc, sfy_loga)
rm(sfy_thimc)

# combine SFY_FALL and SFY (all part of SFY mainstem at CWS site)
sfy <- bind_rows(sfy, sfy_fallck)
rm(sfy_fallck)

# rm remington and grizzly
rm(mfy_remmington) 
rm(mfy_grizzly_ck)

# combine deer:
deer_clearck <- bind_rows(deer_clearck, deer_clec)
rm(deer_clec)

rm(sfy_humbug)

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
sites<-ls(pattern = "deer|nfy|mfy|sfy") # need to make a list of site names

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
sitename <- "yuba"
paste0("put bamlist_mrg_",site,"*",bamNo,"k*")

# cd projects/rana_rapture/MERGED/bamlists/POPS
#paste0("put bamlist_mrg_",site,"_*",bamNo,"k*")

# 7. BASH: PCA CALC SITES ------------------------------------------------

# Use angsd to run pca_calc_sites script: 
lsite<- tolower(site)
paste0("sbatch -p high -t 24:00:00 03_pca_calc_sites.sh bamlist_mrg_",site,"_",bamNo,"k", " ", lsite, "_",bamNo,"k", " bamlists/POPS")

# 8a. BASH: Thetas/SFS ---------------------------------------------------

## in bamlists/subpops:

# ls bamlist_mrg_mfy* | grep "k$" | sed 's/bamlist_mrg_//g' > thetalist_yuba
# ls bamlist_mrg_sfy* | grep "k$" | sed 's/bamlist_mrg_//g' >> thetalist_yuba
#ls bamlist_mrg_nfy* | grep "k$" | sed 's/bamlist_mrg_//g' >> thetalist_yuba
#ls bamlist_mrg_deer* | grep "k$" | sed 's/bamlist_mrg_//g' >> thetalist_yuba
## mv list out:
# mv thetalist_yuba ../../

# sbatch -p high -t 12:00:00 06_theta_sfs.sh thetalist_yuba

# 8. BASH: Calc Pairwise FST ---------------------------------------------

sites<-ls(pattern = "mfy|nfy|sfy|deer") 

# get all unique pairs:
sitepairs <- combn(x = sites, m = 2)

for(i in 1:ncol(sitepairs)){
  cat(paste0("sbatch -p high -t 12:00:00 07_get_fst.sh ",sitepairs[,i][1],"_",bamNo,"k", " ",sitepairs[,i][2],"_",bamNo,"k"))
  cat("\n")
}


# FOR POPS
pops <- c("RABO_bear_100k", "RABO_eel_100k", "RABO_fea_100k", "RABO_mfa_100k",
          "RABO_mfa_hydro_100k","RABO_mfy_100k","RABO_nfa_100k", 
          "RABO_nff_100k", "RABO_nfmfa_100k", "RABO_nfy_100k", "RABO_rub_100k",
          "RABO_sfa_100k", "RABO_sfy_100k", "RABO_trin_100k", "RABO_tuo_100k",
          "RABO_vanduz_100k", "RABO_yuba_100k")
sitepairs <- combn(x=pops, m=2)
sitepairs
for(i in 1:ncol(sitepairs)){
  cat(paste0("sbatch -p high -t 12:00:00 07_get_fst.sh ",sitepairs[,i][1], " ",sitepairs[,i][2]))
  cat("\n")
}

#### POPS
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_eel_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_fea_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_mfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_mfa_hydro_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_mfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_nfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_nff_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_nfmfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_bear_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_fea_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_mfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_mfa_hydro_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_mfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_nfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_nff_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_nfmfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_eel_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_mfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_mfa_hydro_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_mfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_nfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_nfmfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_fea_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_mfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_nfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_nff_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_nfmfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_mfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_nfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_nff_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_nfmfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfa_hydro_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_nfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_nff_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_nfmfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_mfy_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_nff_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_nfmfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfa_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_nfmfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nff_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfmfa_100k RABO_nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfmfa_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfmfa_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfmfa_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfmfa_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfmfa_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfmfa_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfmfa_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfy_100k RABO_rub_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfy_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfy_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfy_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfy_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_nfy_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_rub_100k RABO_sfa_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_rub_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_rub_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_rub_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_rub_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_rub_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_sfa_100k RABO_sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_sfa_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_sfa_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_sfa_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_sfa_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_sfy_100k RABO_trin_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_sfy_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_sfy_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_trin_100k RABO_tuo_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_trin_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_trin_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_tuo_100k RABO_vanduz_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_tuo_100k RABO_yuba_100k
sbatch -p high -t 12:00:00 07_get_fst.sh RABO_vanduz_100k RABO_yuba_100k


# can check datestamp & sort by time: ls -lt FILE*
### SUBPOPS
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k mfy_oregck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k mfy_oregck_flt_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k mfy_us_oh_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k nfy_slate_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k sfy_loga_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k sfy_mcki_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh deer_clearck_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k mfy_us_oh_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k nfy_slate_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k sfy_loga_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k sfy_mcki_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k mfy_us_oh_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k nfy_slate_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k sfy_loga_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k sfy_mcki_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_oregck_flt_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k nfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k nfy_slate_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k sfy_loga_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k sfy_mcki_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh mfy_us_oh_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k nfy_slate_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k sfy_loga_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k sfy_mcki_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_slate_100k sfy_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_slate_100k sfy_loga_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_slate_100k sfy_mcki_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_slate_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_slate_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_slate_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_slate_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh nfy_slate_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_100k sfy_loga_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_100k sfy_mcki_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_loga_100k sfy_mcki_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_loga_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_loga_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_loga_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_loga_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_loga_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_mcki_100k sfy_misc_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_mcki_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_mcki_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_mcki_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_mcki_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_misc_100k sfy_rockck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_misc_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_misc_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_misc_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_rockck_100k sfy_scotchman_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_rockck_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_rockck_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_scotchman_100k sfy_shadyck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_scotchman_100k sfy_springck_100k
sbatch -p high -t 12:00:00 07_get_fst.sh sfy_shadyck_100k sfy_springck_100k

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




