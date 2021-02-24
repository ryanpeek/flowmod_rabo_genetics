# COMBINE THE WAVELET AND COLWELL DATA AND PLOT

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(WaveletComp) # for wavelet analysis
library(hydrostats) # for seasonality colwell analysis
library(viridis)
library(ggforce)
library(ggrepel)
library(sf)
library(ggridges) 
library(ggthemes)


# FUNCTIONS ---------------------------------------------------------------

source("scripts/functions/f_USGS_get_raw_dat.R")
source("scripts/functions/f_doy.R")
source("scripts/functions/f_get_usgs_dv.R")
source("scripts/functions/f_dv_facet_zoom.R")

# function to read in RData and rename on the fly:
loadRData <- function(fileName){
  load(fileName) # loads file
  get(ls()[ls() != "fileName"]) # returns but removes the name
}

# LOAD & MERGE DATA ------------------------------------------------------

# List files for a given site/river for MAINSTEM only (_dv)
(files <- list.files(path = "models/dv_flow/", pattern = "*rda"))

# pre allocate a list to use:
wave.df <- list()
cole.df <- list()

for(i in seq_along(files)){
  print(paste0("Processing: ",files[i]))
  load(paste0("models/dv_flow/",files[i])) 
  wave.df[[i]]<-bind_rows(predict.mf, wave.df)
  cole.df[[i]]<- bind_rows(season.df, cole.df)
  print(paste0("Completed: ", files[i]))
}

# WAVELETS
wave.df <- map_df(wave.df, bind_rows)
# need to add additional files (RUB_PCWA):
wave.df <- bind_rows(wave.df, p.mf)
rm(p.mf)

# COLWELL
cole.df <- map_df(cole.df, bind_rows) %>% distinct(site, .keep_all = T)
# add additional single RUB PCWA
s.mp <- s.mp %>% rename(MP=MP_metric)
cole.df <- bind_rows(cole.df, s.mp)
rm(s.mp, season.df, predict.mf)

# BIND SINGLE FILES -------------------------------------------------------

# fix and bind MFY OC
load("models/MFY_OC_11409400_usgs_wavelet_colwell.rda") # ORCK below Log Cabin 1968-2016

mfy_oc.season <- rename(mfy.season, MP=MP_metric) %>% 
  mutate(site="ORCK_11409400_1968_2016")
cole.df <- bind_rows(cole.df, mfy_oc.season)
rm(mfy_oc.season, mfy.season)

plotMFY_oc <- plotMFY %>% 
  mutate(site="ORCK_11409400_1968_2016")
wave.df <- bind_rows(wave.df, plotMFY_oc)
rm(plotMFY_oc, plotMFY)

# fix and bind MFY OC US
load("models/MFY_OC_US_11409300_usgs_wavelet_colwell.rda") # ORCK upstream Log Cabin 1967-2000

mfy_oc_us.season <- rename(mfy_oc_us.season, MP=MP_metric) %>% 
  mutate(site="ORCK_11409300_1967_2000")
cole.df <- bind_rows(cole.df, mfy_oc_us.season)
rm(mfy_oc_us.season)

plotMFY_oc_us <- plotMFY_oc_us %>% 
  mutate(site="ORCK_11409300_1967_2000")
wave.df <- bind_rows(wave.df, plotMFY_oc_us)
rm(plotMFY_oc_us)



# DISTILL AND SAVE --------------------------------------------------------

# reduce to max period and max power value
(wavelet <- wave.df %>% group_by(site) %>% 
   filter(round(Period,0)==12) %>%  # filter to period=12
   distinct() %>% # get rid of duplicates
   summarize(
     avgPower = mean(Power.avg, na.rm = T),
     maxPowerAvg=max(Power.avg, na.rm = T),
     maxPeriod=Period[which.max(Power.avg)]))

combined_df <- left_join(wavelet, cole.df, by="site") %>% 
  mutate(river=str_extract(site, "[A-Z]+"),
         gageID=str_extract(site, "[:digit:]{8}"),
         year_sta=str_extract(site, "(?<=[:digit:]{8}_)[:digit:]{4}"),
         year_end=str_extract(site, "[:digit:]{4}$"))

# fix RUB
combined_df$year_sta[30]<-2009

# remove rows with fewer than 10 yrs or duplicate gage data (SFY, MFY, RUB PCWA):
combined_df<- combined_df[-c(17, 18, 30, 40),]

# add reg col
combined_df$reg<-ifelse(combined_df$MP>0.74, "N", "Y")  

# SAVE
#save(combined_df, file = "models/wavelet_colwell_combined.rda")
#write_csv(combined_df, path = "models/wavelet_colwell_combined.csv")
