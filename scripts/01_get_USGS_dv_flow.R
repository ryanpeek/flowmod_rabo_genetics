# GET CDEC & USGS flow

# NOTES:

# Use the file in: "data/streamgages_06_usgs_update.kmz" to view USGS stations, then pull the Gage ID for each watershed. 

# Plot data and look for eveness and gaps, filter and export. Some are pre-dam data, so designate as such.

# use those datasets in next step to calculate seasonality and colwell's MP using wavelet and hydrostats libraries.


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(WaveletComp) # for wavelet analysis
library(hydrostats) # for seasonality colwell analysis
library(viridis)
library(ggforce)
library(ggrepel)
library(sf)
#library(plotly)

# FUNCTIONS ---------------------------------------------------------------

source("scripts/functions/f_getCDEC.R")
# if script doesn't work use historical CSV download:
# "http://cdec.water.ca.gov/cgi-progs/queryCSV"

source("scripts/functions/f_USGS_get_raw_dat.R")
source("scripts/functions/f_doy.R")
source("scripts/functions/f_get_usgs_dv.R")
source("scripts/functions/f_dv_facet_zoom.R")

# TEST FUNCTION ---------------------------------------------------------

get_usgs_dv(gage = 11422500, site = "BEAR", facetYr = 1975, filterYrs = T, savedat = T, saveplot = F)

## FUCKING AMAZINGW!!!

# READ KMZ/KMLS -----------------------------------------------------------
# kmz works on linux only (at the moment)
#if(Sys.info()[1] == "Linux")
#  gages_kmz <- st_read("data/streamgages_06_usgs_update.kmz")

# otherwise use "unzip" option:
gages_kml <- st_read(unzip("data/streamgages_06_usgs_update.kmz"))

# remove the zipped version (usually just called "doc.kml")
file.remove(list.files(pattern = "*.kml",recursive = F))

# plot(st_geometry(gages_kml))

# READ HUCS ----------------------------------------------------------------

hucs <- st_read(unzip("data/shps/HUC8_named_westcoast.zip"), quiet = F) %>%
  st_transform(crs=4326) %>%
  mutate(lon=map_dbl(geometry, ~st_centroid(.x)[[1]]),
         lat=map_dbl(geometry, ~st_centroid(.x)[[2]])) # add centroid values

file.remove(list.files(pattern = "HUC8_named_westcoast*",recursive = F))

#huc_tst<-filter(hucs, grepl("Bear", x = HU_8_NAME)) # find out HUC number

# INTERSECT GAGES WITH HUCS -----------------------------------------------

# intersect a watershed to view all gages within that watershed:
huc_s <- filter(hucs, grepl("South Fork American", x=HU_8_NAME))
plot(st_geometry(huc_s))

# intersect (clip)
gages_s<-st_intersection(gages_kml, huc_s)

# plot with sf
plot(st_geometry(huc_s), border = "blue")
plot(st_geometry(gages_s), add=T)

# plot with ggplot
ggplot() + 
  geom_sf(data=huc_s, color="blue", fill=NA) +
  geom_sf(data=gages_s, color="orange", pch=16)

# simple leaflet:
library(leaflet)
leaflet() %>% 
  addTiles() %>% 
  addPolygons(data=huc_s, fillOpacity = 0.1, fill = NA) %>% 
  addCircleMarkers(data = gages_s, color = "orange", fillColor = "orange", fillOpacity = 0.8, radius = 6, stroke=T, weight = 2)

# USGS: dataRetrieval ------------------------------------------------------

# # check what data exists: can use siteNumber, stateCd, huc, bBox, or countyCd.
# library(dataRetrieval)
# 
# data_ca <- whatNWISdata(huc=18020126)
# nrow(data_ca)
# 
# head(data_ca)
# 
# # unique:
# sites_ca <- unique(data_ca$site_no)
# length(sites_ca)
# 
# # updated streams (temp: 00010, Q: 00065)
# #data_ca_str <- whatNWISdata(huc=18020126, siteType="ST", parameterCd="00065")
# 
# 
# library(maps)
# map('state', regions='california')
# points(x=data_ca$dec_long_va, 
#        y=data_ca$dec_lat_va)
# 

# NFA ----------------------------------------------------------------

# USGS: NFA

# at Shirttail
# 11426500 (1911-1941)
get_usgs_dv(gage = 11426500, site = "NFA", facetYr = 1934, filterYrs = F, savedat = T, saveplot = F)

# at NF Clementine
# 11427000 (1941-2017)
get_usgs_dv(gage = 11427000, site = "NFA", facetYr = 1974, filterYrs = F, savedat = T, saveplot = F)

# quick wavelet analysis:
load("data/flows/NFA_dv_USGS_11426500_1911_1941.rda")
nfa <- filename

site <- "NFA"

nfa.w <- analyze.wavelet(nfa, my.series = 4, dt = 1/30) # usgs flow is in col 4

wt.image(nfa.w, label.time.axis = T, show.date=TRUE, 
         date.format = "%Y",
         main = paste0(site, ": Seasonality of Daily Flow"),
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (years)", periodlab = "period (months)")

wt.avg(nfa.w, show.siglvl = T, siglvl = c(.001, 0.01, 0.05))

plotNFA<-nfa.w[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame

ggplot() + 
  geom_line(data=plotNFA, aes(x=Period, y=Power.avg)) +
  geom_point(data=plotNFA, aes(x=Period,y=Power.avg), 
             col=ifelse(plotNFA$Power.avg.pval<0.05, "red", "blue")) + theme_bw() +
  scale_x_continuous(breaks=seq(0,64, 2), limits = c(0,64))

# quick colwell analysis of seasonality:
## From Tonkin et al 2017: M (Contingency) as metric of seasonality
## To standardize the role of seasonality in relation to overall predictability,
## we divided (M) by overall predictability (the sum of (M) and constancy (C)

# standardize
df.nfa <- nfa %>%  
  rename(Date=date, Q=flow_cfs)

# analyze
nfa.Col <- hydrostats::Colwells(df.nfa)
(nfa.season <- tibble(site=c(site), MP_metric=c(nfa.Col$MP)))

# save models
save(plotNFA, nfa.Col, nfa.season,  file = paste0("models/",site,"_usgs_wavelet_colwell.rda"))

# IV DATA

# load IV data
load("data/flows/NFA_iv_USGS_updated_2017-02-22.rda")
NFA_iv <- add_WYD(NFA_iv, datecolumn = "datetime")

# plot
ggplot() + geom_line(data=NFA_iv, aes(x=datetime, y=flow_cfs), show.legend = F, color="darkblue") + ylab("log[Discharge] (cms)") +xlab("")+ scale_y_log10()+
  facet_zoom(x = WY == 2012, shrink=T)# + ylim(c(0,500))

#ggsave(filename = "figs/NFA_log_flow_facet_zoom.png",width = 9, height = 7, units = "in")



# MFA ---------------------------------------------------------------------

# At Quarry
# 11433500 (1911-1986) # dam operational starting in 1965/66? (1911-1964, 1965-1986)
get_usgs_dv(gage = 11433500, site = "MFA", facetYr = 1965, filterYrs = T, savedat = T, saveplot = T)

dv_zoom(dat_dv, "MFA 11433500", scalelog = T, zoomYr = 1969)

# Main Bar Canyon Ck (Near American Canyon)
# 11433420 (1972-1986)
get_usgs_dv(gage = 11433420, site = "MFA_MBC", facetYr = 1976, filterYrs = F, savedat = T, saveplot = F)

# Canyon C NR Georgetown
# 11433400 (1966-1979) # dam operational starting in 1966?
get_usgs_dv(gage = 11433400, site = "MFA_CCk", facetYr = 1969, filterYrs = F, savedat = T, saveplot = F)

# MFA at Tunnel Chute
# 11433300 (1958-2016) # dam in 1965 so trim to >1964
get_usgs_dv(gage = 11433300, site = "MFA", facetYr = 1966, filterYrs = T, savedat = T, saveplot = F)

dat_dv <- dat_dv[-1,]
dat_dv$flow_cfs <- as.numeric(dat_dv$flow_cfs)
dv_zoom(dat_dv, "MFA 11433300", scalelog = T, zoomYr = 1965)

# MFA at Oxbow PH
# 11433212 (1973-2005)
get_usgs_dv(gage = 11433212, site = "MFA", facetYr = 1990, filterYrs = F, savedat = T, saveplot = F)

# MFA at Ralston PH
# 11427765 (1973-2005)
get_usgs_dv(gage = 11427765, site = "MFA", facetYr = 1990, filterYrs = F, savedat = T, saveplot = F)

# MFA at INTERBAY (bypass)
# 11427770 (1965-2016)
get_usgs_dv(gage = 11427770, site = "MFA", facetYr = 1990, filterYrs = F, savedat = T, saveplot = F)

# quick wavelet analysis:
mfa.w <- analyze.wavelet(MFA_dv, my.series = 3, dt = 1/30)

wt.image(mfa.w, main = "MFA Seasonality of Daily Flow",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (days)", periodlab = "period (months)")

wt.avg(mfa.w)

plotMFA<-mfa.w[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame

ggplot() + geom_line(data=plotMFA, aes(x=Period, y=Power.avg))+geom_point(data=plotMFA, aes(x=Period,y=Power.avg), col=ifelse(plotMFA$Power.avg.pval<0.05, "red", "blue")) + scale_x_continuous(breaks=seq(0,64, 2), limits = c(0,64))

# quick colwell analysis of seasonality:
## From Tonkin et al 2017: M (Contingency) as metric of seasonality
## To standardize the role of seasonality in relation to overall predictability,
## we divided (M) by overall predictability (the sum of (M) and constancy (C)

# standardize
df.mfa <- MFA_dv %>%  
  rename(Date=date, Q=flow_avg_cfs) %>% as.data.frame()

# analyze
Col.mfa<-hydrostats::Colwells(df.mfa)
(seasonality <- tibble(site=c("mfa"), MP_metric=c(Col.mfa$MP)))

# RUBICON -----------------------------------------------------------------
# Long Canyon
# 11433100 (1960-1992)
get_usgs_dv(gage = 11433100, site = "RUB_LC", facetYr = 1970, filterYrs = T, savedat = T, saveplot = F)

# At PH
# 11433200 (1958-84) POST/PRE DAM (64/65)
get_usgs_dv(gage = 11433200, site = "RUB", facetYr = 1965, filterYrs = T, savedat = T, saveplot = F)

dv_zoom(dat_dv, "MFA 11433200", scalelog = T, zoomYr = 1965)

# US at Crossing # PRE_dam in 1966
# 11431000 (1943-1964)
get_usgs_dv(gage = 11431000, site = "RUB", facetYr = 1964, filterYrs = T, savedat = F, saveplot = F)

# wavelet analysis
rub.w <- analyze.wavelet(RUB_dv, my.series = 3, dt = 1/30)

wt.image(rub.w, main = "RUB Seasonality of Daily Flow",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (days)", periodlab = "period (months)")

wt.avg(rub.w)

plotRUB<-rub.w[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame

ggplot() + geom_line(data=plotRUB, aes(x=Period, y=Power.avg))+geom_point(data=plotRUB, aes(x=Period,y=Power.avg), col=ifelse(plotRUB$Power.avg.pval<0.05, "red", "blue")) + scale_x_continuous(breaks=seq(0,64, 2), limits = c(0,64))

# quick colwell analysis of seasonality:
## From Tonkin et al 2017: M (Contingency) as metric of seasonality
## To standardize the role of seasonality in relation to overall predictability,
## we divided (M) by overall predictability (the sum of (M) and constancy (C)

# standardize
df.RUB <- RUB_dv %>%  
  rename(Date=date, Q=flow_avg_cfs) %>% as.data.frame()

# analyze
Col.RUB<-hydrostats::Colwells(df.RUB)
(seasonality <- tibble(site=c("RUB"), MP_metric=c(Col.RUB$MP)))

# NFMFA -------------------------------------------------------------------

# USGS:
# 11433260 (1965-85)
get_usgs_dv(gage = 11433260, site = "NFMFA", facetYr = 1970, filterYrs = F, savedat = T, saveplot = F)

# SFA -------------------------------------------------------------------
# Chili Bar in 1964, Slab 1967

# At Lotus
# 11445500 (1951-95) # flow reg pattern appears 1961
get_usgs_dv(gage = 11445500, site = "SFA", facetYr = 1962, filterYrs = T, savedat = T, saveplot = F)

dv_zoom(dat_dv, "SFA", scalelog = T, zoomYr = 1961)

# At Coloma (PRE REG)
# 11445000 (1929-1941)
get_usgs_dv(gage = 11445000, site = "SFA", facetYr = 1935, filterYrs = F, savedat = T, saveplot = F)

# Near Placerville
# 11444500 (1911-2005) 1911-1920, 1965-2016
get_usgs_dv(gage = 11444500, site = "SFA", facetYr = 1915, filterYrs = T, savedat = T, saveplot = F)

# SFA Below Camino Dam
# 11443500 (1922-2016) # dam in 1961, flows alter/baseflow begin 1959
get_usgs_dv(gage = 11443500, site = "SFA", facetYr = 1958, filterYrs = T, savedat = T, saveplot = F)

dv_zoom(dat_dv, "SFA", scalelog = T, zoomYr = 1961)


# Silver Ck (Camino) Near Conf
# 11442000 (1922-1961)
get_usgs_dv(gage = 11442000, site = "SFA_SILV", facetYr = 1960, filterYrs = T, savedat = T, saveplot = F)

# Silver Ck DS of Camino
# 11441900 (1960-2016)
get_usgs_dv(gage = 11441900, site = "SFA_SILV", facetYr = 1960, filterYrs = T, savedat = T, saveplot = F)


# BEAR --------------------------------------------------------------------

# USGS GAGE INFO:

# Below Rollins:
# 11422500 (1965 to 2016) # this one best [filter to after 1965]
get_usgs_dv(gage = 11422500, site = "BEAR", facetYr = 1975, filterYrs = T, savedat = T, saveplot = F)

## 11421900 (1986 to 2016) # better dataset, shows hydro/bypass fx
get_usgs_dv(gage = 11421900, site = "BEAR", facetYr = 1995, filterYrs = T, savedat = T, saveplot = F)

### 11422000 (1912 to 2016, missing 1953-1964), not ideal

# Below 49:
## 11423000 (1941 to 1967) # good for early period
get_usgs_dv(gage = 11423000, site = "BEAR", facetYr = 1945, filterYrs = F, savedat = T, saveplot = F)

# Above Rollins at Dutch Flat 1:
# 11421790 (1965 to 2016) # DUTCH AFTERBAY: good but very low/zero flow
get_usgs_dv(gage = 11421790, site = "BEAR", facetYr = 1995, filterYrs = F, savedat = T, saveplot = F)

# Blue Canyon at Drum PH:
# 11421770 (1966 to 2016) # DRUM AFTERBAY: good dataset, shows min flows after 2000
get_usgs_dv(gage = 11421770, site = "BEAR", facetYr = 1995, filterYrs = F, savedat = T, saveplot = F)

## 11414196 (1981 to 2013) # DRUM PH 2: combined flow drum 1 and 2 pp
get_usgs_dv(gage = 11414196, site = "BEAR", facetYr = 1995, filterYrs = F, savedat = T, saveplot = F)



# quick wavelet analysis:
bear.w <- analyze.wavelet(bear, my.series = 4, dt = 1/30) # usgs flow is in col 4

wt.image(bear.w, main = paste0(site, ": Seasonality of Daily Flow"),
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (days)", periodlab = "period (months)")

wt.avg(bear.w)

plotBEAR<-bear.w[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame

ggplot() + geom_line(data=plotBEAR, aes(x=Period, y=Power.avg))+geom_point(data=plotBEAR, aes(x=Period,y=Power.avg), col=ifelse(plotBEAR$Power.avg.pval<0.05, "red", "blue")) + scale_x_continuous(breaks=seq(0,64, 2), limits = c(0,64))

# quick colwell analysis of seasonality:
## From Tonkin et al 2017: M (Contingency) as metric of seasonality
## To standardize the role of seasonality in relation to overall predictability,
## we divided (M) by overall predictability (the sum of (M) and constancy (C)

# standardize
df.bear <- bear %>%  
  rename(Date=date, Q=flow_cfs)

# analyze
bear.Col <- hydrostats::Colwells(df.bear)
(bear.season <- tibble(site=c(site), MP_metric=c(bear.Col$MP)))

# save models
save(plotBEAR, bear.Col, bear.season,  file = paste0("models/",site,"_usgs_wavelet_colwell.rda"))

# NFY ---------------------------------------------------------------

# NFY
# 11410500 (1923-1944) at Sierra City
get_usgs_dv(gage = 11410500, site = "NFY", facetYr = 1934, filterYrs = F, savedat = T, saveplot = F)
# 11413000 (1930-2017) # goodyears
get_usgs_dv(gage = 11413000, site = "NFY", facetYr = 1964, filterYrs = F, savedat = T, saveplot = F)
# 11413100 (1968-1987) # NFY above slate
get_usgs_dv(gage = 11413100, site = "NFY", facetYr = 1974, filterYrs = F, savedat = T, saveplot = F)

# Slate Ck
# 11413300 (1960-2005) below diver dam
get_usgs_dv(gage = 11413300, site = "NFY_SC", facetYr = 1974, filterYrs = F, savedat = T, saveplot = F)

# quick wavelet analysis:

nfy.w <- analyze.wavelet(nfy, my.series = 4, dt = 1/30) # usgs flow is in col 4

wt.image(nfy.w, main = paste0(site, ": Seasonality of Daily Flow"),
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (days)", periodlab = "period (months)")

wt.avg(nfy.w)

plotNFY<-nfy.w[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame

summary(plotNFY)

ggplot() + geom_line(data=plotNFY, aes(x=Period, y=Power.avg)) + 
  geom_point(data=plotNFY, aes(x=Period,y=Power.avg), 
             col=ifelse(plotNFY$Power.avg.pval<0.05, "red", "blue")) + 
  scale_x_continuous(breaks=seq(0,36, 2), limits = c(0,36)) + 
  theme_bw()

# quick colwell analysis of seasonality:
## From Tonkin et al 2017: M (Contingency) as metric of seasonality
## To standardize the role of seasonality in relation to overall predictability,
## we divided (M) by overall predictability (the sum of (M) and constancy (C)

# standardize
df.nfy <- nfy %>%  
  rename(Date=date, Q=flow_cfs)

# analyze
nfy.Col <- hydrostats::Colwells(df.nfy)
(nfy.season <- tibble(site=c(site), MP_metric=c(nfy.Col$MP)))

# save models
save(plotNFY, nfy.Col, nfy.season,  file = paste0("models/",site,"_",gage,"_usgs_wavelet_colwell.rda"))

# MFY ---------------------------------------------------------------

# USGS: MFY

# 11409400 (1968-2016) # at Oregon Ck below Log Cabin Dam (oc)
get_usgs_dv(gage = 11409400, site = "MFY_OC", facetYr = 1984, filterYrs = F, savedat = T, saveplot = F)

# load("/Users/ryanpeek/Documents/github/rabo_regulation/data/flows/MFY_OC_dv_USGS_11409400_1968_2016.rda")
# dat_dv <- dat_dv[-1,]
# dat_dv$flow_cfs <- as.numeric(dat_dv$flow_cfs)
# dv_zoom(dat_dv, "MFY_OC", scalelog = T, zoomYr = 1990)
# save(dat_dv, file = "data/flows/MFY_OC_dv_USGS_11409400_1968_2016.rda")

# 11409350 (1988-2016) # check (oc2) # tunnel? # drops to zero
#get_usgs_dv(gage = 11409350, site = "MFY_OC", facetYr = 1994, filterYrs = F, savedat = T, saveplot = F)
# dat_dv <- dat_dv[-1,]
# dat_dv$flow_cfs <- as.numeric(dat_dv$flow_cfs)
# dv_zoom(dat_dv, "MFY_OC", scalelog = T, zoomYr = 1995)

# 11409300 (1967-2000) # US Log Cabin (oc_us)
get_usgs_dv(gage = 11409300, site = "MFY_OC_US", facetYr = 1984, filterYrs = F, savedat = T, saveplot = F)

# MFY
# 11410000 (1913-1940 & 2001-2005) # ds ORCk at 49 bridge (OH DAM BUILT ~1968, Bowman 1927)
get_usgs_dv(gage = 11410000, site = "MFY", facetYr = 1934, filterYrs = T, savedat = T, saveplot = F)
get_usgs_dv(gage = 11410000, site = "MFY", facetYr = 2004, filterYrs = T, savedat = T, saveplot = F)

# 11409000 (1941-1969)
get_usgs_dv(gage = 11409000, site = "MFY", facetYr = 1964, filterYrs = F, savedat = T, saveplot = F)

# 11408880 (1968-2016)
get_usgs_dv(gage = 11408880, site = "MFY", facetYr = 1974, filterYrs = F, savedat = T, saveplot = F)

# 11408550 (1987-2016) # below Milton (built 1928)
get_usgs_dv(gage = 11408550, site = "MFY", facetYr = 1994, filterYrs = F, savedat = T, saveplot = F)


# quick wavelet analysis:
mfy<-filename
site <- "MFY_OC_US"
gage <- "11409300"
mfy.w <- analyze.wavelet(mfy, my.series = 4, dt = 1/30) # usgs flow is in col 4

wt.image(mfy.w, main = paste0(site, ": Seasonality of Daily Flow"),
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (days)", periodlab = "period (months)")

wt.avg(mfy.w)

plotMFY<-mfy.w[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame

summary(plotMFY)

ggplot() + geom_line(data=plotMFY, aes(x=Period, y=Power.avg)) + 
  geom_point(data=plotMFY, aes(x=Period,y=Power.avg), 
             col=ifelse(plotMFY$Power.avg.pval<0.05, "red", "blue")) + 
  scale_x_continuous(breaks=seq(0,36, 2), limits = c(0,36)) + 
  theme_bw()

# quick colwell analysis of seasonality:
## From Tonkin et al 2017: M (Contingency) as metric of seasonality
## To standardize the role of seasonality in relation to overall predictability,
## we divided (M) by overall predictability (the sum of (M) and constancy (C)

# standardize
df.mfy <- mfy %>%  
  rename(Date=date, Q=flow_cfs)

# analyze
mfy.Col <- hydrostats::Colwells(df.mfy)
(mfy.season <- tibble(site=c(site), MP_metric=c(mfy.Col$MP)))

# rename
plotMFY_oc_us <- plotMFY
mfy_oc_us.season <- mfy.season

# save models
save(plotMFY_oc_us, mfy_oc_us.season,  file = paste0("models/",site,"_",gage,"_usgs_wavelet_colwell.rda"))

# save models
save(plotMFY, mfy.season,  file = paste0("models/",site,"_",gage,"_usgs_wavelet_colwell.rda"))

# SFY ---------------------------------------------------------------

# Langs Crossing # spaulding in early 1913
# 11414250 (1965-2016)
get_usgs_dv(gage = 11414250, site = "SFY", facetYr = 1994, filterYrs = F, savedat = T, saveplot = F)

dv_zoom(dat_dv, "MFA 11433200", scalelog = T, zoomYr = 1965)

# Canyon Ck below Bowman
# 11416500 (1927-2016) not great because very little flow
#get_usgs_dv(gage = 11416500, site = "SFY", facetYr = 1994, filterYrs = F, savedat = T, saveplot = F)

# S Yuba at Washington
# 11417000 (1942-1972)  # gap between '54-56
get_usgs_dv(gage = 11417000, site = "SFY", facetYr = 1964, filterYrs = T, savedat = T, saveplot = F)

# Poorman Ck
# 11417100 (1942-1972)

# SFY at Jones Bar
# 11417500 (1940-2016) # missing 1949 to 1958
get_usgs_dv(gage = 11417500, site = "SFY", facetYr = 1974, filterYrs = T, savedat = T, saveplot = F)


# quick wavelet analysis:

sfy.w <- analyze.wavelet(sfy, my.series = 4, dt = 1/30) # usgs flow is in col 4

wt.image(sfy.w, main = paste0(site, ": Seasonality of Daily Flow"),
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (days)", periodlab = "period (months)")

wt.avg(sfy.w)

plotSFY<-sfy.w[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame

summary(plotSFY)

ggplot() + geom_line(data=plotSFY, aes(x=Period, y=Power.avg)) + 
  geom_point(data=plotSFY, aes(x=Period,y=Power.avg), 
             col=ifelse(plotSFY$Power.avg.pval<0.05, "red", "blue")) + 
  scale_x_continuous(breaks=seq(0,36, 2), limits = c(0,36)) + 
  theme_bw()

# quick colwell analysis of seasonality:
## From Tonkin et al 2017: M (Contingency) as metric of seasonality
## To standardize the role of seasonality in relation to overall predictability,
## we divided (M) by overall predictability (the sum of (M) and constancy (C)

# standardize
df.sfy <- sfy %>%  
  rename(Date=date, Q=flow_cfs)

# analyze
sfy.Col <- hydrostats::Colwells(df.sfy)
(sfy.season <- tibble(site=c(site), MP_metric=c(sfy.Col$MP)))

# save models
save(plotSFY, sfy.Col, sfy.season,  file = paste0("models/",site,"_",gage,"_usgs_wavelet_colwell.rda"))

