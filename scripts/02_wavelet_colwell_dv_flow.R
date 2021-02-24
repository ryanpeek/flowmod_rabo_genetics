# GET CDEC & USGS flow

# NOTES:

# Run the data_01_USGS_dv_flow.R script to pull data

# Merge the data into one giant database and add info for each site, reg/unreg, pre-dam, post-dam

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(WaveletComp) # for wavelet analysis
library(hydrostats) # for seasonality colwell analysis
library(viridis)
library(ggforce)
library(ggrepel)
library(sf)
library(ggridges) # joyplots
#library(plotly)

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

#mf <- loadRData(paste0("data/flows/",files[1]))
#get_usgs_dv(gage = 11422500, site = "BEAR", facetYr = 1975, filterYrs = T, savedat = T, saveplot = F)

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

# # simple leaflet:
# library(leaflet)
# leaflet() %>% 
#   addTiles() %>% 
#   addPolygons(data=huc_s, fillOpacity = 0.1, fill = NA) %>% 
#   addCircleMarkers(data = gages_s, color = "orange", fillColor = "orange", fillOpacity = 0.8, radius = 6, stroke=T, weight = 2)


# LOAD FLOW DATA: SINGLE ------------------------------------------------

river <- "MFA_dv_USGS"

# List files for a given site/river for MAINSTEM only (_dv)
(files <- list.files(path = "data/flows/", pattern = paste0("^",river)))

gageID <- files %>%
  map(stringr::str_split, pattern = "_") %>%
  map_chr(c(1,4))

rivID <- files %>%
  map(stringr::str_split, pattern = "_") %>%
  map_chr(c(1,1))

df <- tibble("fileext"=files,
             "site"=tools::file_path_sans_ext(gsub(files,pattern = "_dv_USGS", replacement = "")),
             rivID, gageID)

# i<-1
# mf <- loadRData(paste0("data/flows/",df[i,1]))


# LOAD FLOW DATA: PURRR ------------------------------------------------------

river <- "SFY_dv_USGS"

# List files for a given site/river for MAINSTEM only (_dv)
(files <- list.files(path = "data/flows/", pattern = paste0("^",river)))

gageID <- files %>%
  map(stringr::str_split, pattern = "_") %>%
  map_chr(c(1,4))

rivID <- files %>%
  map(stringr::str_split, pattern = "_") %>%
  map_chr(c(1,1))

df <- tibble("fileext"=files,
             "site"=tools::file_path_sans_ext(gsub(files,pattern = "_dv_USGS", replacement = "")),
             rivID, gageID)

mf <- data_frame(filename = files,
                   site = tools::file_path_sans_ext(gsub(files,pattern = "_dv_USGS", replacement = ""))) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ loadRData(file.path("data/flows/", .))),
         file_contents = map(file_contents, ~mutate_at(.x, "gageNo", as.character)))# %>% unnest


mf_flat <- data_frame(filename = files,
                 site = tools::file_path_sans_ext(gsub(files,pattern = "_dv_USGS", replacement = ""))) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ loadRData(file.path("data/flows/", .))),
         file_contents = map(file_contents, ~mutate_at(.x, "gageNo", as.character))) %>% unnest

# PREVIEW PLOT ------------------------------------------------------------

tst1 <- filter(mf_flat, gageNo=="11417500")
rivname <- "SFY"

# ridgeplot
ggplot(data=tst1, aes(x=DOWY, y=WY, height = log(flow_cfs), group = WY)) +
  geom_density_ridges(stat = "identity", scale = 2.5, fill="skyblue")

# line plot w intercept of Ralston Dam install
ggplot() + 
  geom_line(data=mf_flat, aes(x = date, y = flow_cfs, group=site, color=flow_cfs)) + 
  facet_grid(site~.) + scale_y_log10() + 
  geom_vline(xintercept = ymd("1966-01-01"), lty=2, col="maroon") +
  theme_bw() + scale_color_viridis()

# SIMPLE FACET_ZOOM
source("scripts/functions/f_dv_facet_zoom.R")
dv_zoom(tst1, rivname, scalelog = T, zoomYr = 1990)

# FULL FACET ZOOM

facet_yr <- 1990

# facet zoom
ggplot() + 
  geom_line(data=tst1, aes(x=date, y=flow_cfs, color=DOWY)) +
  scale_y_log10(labels = scales::comma) + 
  #ylim(0,10000)+
  theme_bw(base_family = "Roboto Condensed")+
  ylab("Discharge (cfs)") + xlab("") +
  ggtitle(paste0(rivname,": USGS ", 
                 tst1$gageNo[1], " (", length(unique((tst1$WY)))," years)")) +
  scale_color_viridis_c(option = "D") +
  ggforce::facet_zoom(x = WY == facet_yr, shrink=T)

#ggsave(filename = paste0("figs/", tst1$site[1],"_facet",facet_yr,".png"), width = 9, height = 7, units = "in", dpi = 200)

# WAVELET LOOP ------------------------------------------------------------

# pre allocate a list to use:
wt.list <- list()
for(i in seq_along(mf[[1]])){
  print(paste0("Processing: ",mf$filename[i]))
  #wt.list[[i]] <- tibble(file=mf$filename[i], site=mf$site[i]) # for testing
  wt.list[[i]] <- analyze.wavelet(mf[[3]][[i]], my.series = 4, dt = 1/30)
  print(paste0("Completed: ", mf$filename[i]))
}

# set names in list
wt.list <- wt.list %>% set_names(df$site)

for(i in seq_along(wt.list)){
  
  pdf(file = paste0("figs/",mf[[2]][i],".pdf"))
  
  wt.image(wt.list[[i]], main = paste0(names(wt.list[i]),
                                       ": Seasonality of Daily Flow"),
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (days)", periodlab = "period (months)")
  
  dev.off()
}

# pre allocate data_frame
predict.mf <- list()

for(i in seq_along(wt.list)){
  print(paste0("Processing: ",names(wt.list[i])))
  predict.mf[[i]]<-wt.list[[i]][c("Power.avg","Period","Power.avg.pval")] %>% 
    as.data.frame %>% 
    mutate(site=df$site[i])
}

# flatten all dataframes into one df
predict.mf <- map_df(predict.mf, bind_rows)


# WAVELET SINGLE ----------------------------------------------------------

#wt.mf <- analyze.wavelet(mf[[3]][[3]], my.series = 4, dt = 1/30)

# wt.image(wt.mf, main = paste0(mf$site[3], ": Seasonality of Daily Flow"),
#          legend.params = list(lab = "cross-wavelet power levels"),
#          timelab = "Time (days)", periodlab = "period (months)")

wt.mf <- analyze.wavelet(RUB_dv, my.series = 3, dt = 1/30)
pdf(file = paste0("figs/","RUB_PCWA_2009_2016",".pdf"))
wt.image(wt.mf, main = paste0("RUB_PCWA_2009_2016", ": Seasonality of Daily Flow"),
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time (days)", periodlab = "period (months)")
dev.off()

#wt.avg(wt.mf)
p.mf<-wt.mf[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame %>%
  mutate(site=df$site[3])

p.mf<-wt.mf[c("Power.avg","Period","Power.avg.pval")] %>% as.data.frame %>%
  mutate(site="RUB_PCWA_2009_2016")
 
# ggplot() + geom_line(data=p.mf, aes(x=Period, y=Power.avg)) +
#   geom_point(data=p.mf, aes(x=Period,y=Power.avg), col=ifelse(p.mf$Power.avg.pval<0.05, "red", "blue")) +
#   geom_vline(data=p.mf, xintercept=p.mf$Period[which.max(p.mf$Power.avg)], col="red", lty=2) +
#   annotate(geom = "text", x=4, y=22,
#            label=paste0("Max Power: \n",
#                         round(max(p.mf$Power.avg),digits = 2)))+
#   scale_x_continuous(breaks=seq(0,64, 2), limits = c(0,64)) +
#   theme_bw() +
#   labs(title=paste0(mf$site[3],": Daily Flow"), x="Period (months)")


# SEASONALITY LOOP --------------------------------------------------------

## From Tonkin et al 2017: M (Contingency) as metric of seasonality
## To standardize the role of seasonality in relation to overall predictability,
## we divided (M) by overall predictability (the sum of (M) and constancy (C)
# 
# s.df <- mf[[3]][3] %>% as.data.frame %>%
  # rename(Date=date, Q=flow_cfs)
# s.df <- RUB_dv %>% as.data.frame %>%
  # rename(Date=date, Q=flow_avg_cfs)

# analyze
# s.col <- hydrostats::Colwells(s.df)
# (s.mp <- tibble(site=c("RUB_PCWA_2009_2016"), MP_metric=c(s.col$MP)))

# map over the list of data and calc colwell info
season.list <- data_frame(
  site=mf$site,
  MPdat = map(mf[[3]], rename, Date=date, Q=flow_cfs) %>% 
  map(., hydrostats::Colwells)) #%>% 

season.df <- list()
for(i in seq_along(season.list[[2]])){
  print(paste0("Processing: ",season.list$site[i]))
  season.df[[i]] <- tibble(site=season.list$site[i],
                        MP_metric=season.list[[2]][[i]][7]) %>% flatten
  print(paste0("Completed: ", season.list$site[i]))
}

(season.df <- season.df %>% bind_rows) # flatten it out into df

# SAVE --------------------------------------------------------------------

# see the results:
predict.mf %>% group_by(site) %>% 
  summarize(maxPower=max(Power.avg, na.rm = T),
            maxPeriod=Period[which.max(Power.avg)])

season.df

# save model outputs
save(predict.mf, season.df,  file = paste0("models/", river,"_wavelet_colwell.rda"))

# save singles:
# save(p.mf, s.mp, file = paste0("models/RUB_dv_PCWA_wavelet_colwell.rda"))

