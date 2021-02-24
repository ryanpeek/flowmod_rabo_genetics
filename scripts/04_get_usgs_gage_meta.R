# GET USGS GAGE DATA


# LIBRARIES ---------------------------------------------------------------

library(dataRetrieval)
library(sf)
library(tidyverse)
library(mapview)
library(purrr)
library(lubridate)
library(tidylog)


# GET SITES ---------------------------------------------------------------



# load data (from excel table)
#gage_list <- c(11413000, 11408880, 11409000, 11410000, 11408550, 11414250, 11417500, 11417000, 11426500, 11427000, 11431000, 11433200, 11427765, 11433300, 11427770, 11433500, 11421770, 11423000)

load("models/wavelet_colwell_combined.rda")

combined_df <- combined_df %>% filter(!river %in% c("SFA", "NFMFA", "ORCK")) %>% 
  mutate(site = if_else(site=="MFA_11427765_1973_2016", "RUB_11427765_1973_2016", site),
         river = if_else(site=="RUB_11427765_1973_2016", "RUB", river)) %>% 
  mutate(regtype=case_when(
    river %in% c("RUB","SFY", "MFY", "BEAR") ~ "Bypass",
    river %in% c("NFY", "NFA") ~ "Unregulated",
    river %in% c("MFA") ~ "Hydropeaking"
  )) 

combined_df$regtype <- as.factor(combined_df$regtype)



## DROP SOME SITES TO CLEAN UP PICTURE: all
#gages_to_drop <- c(11422500, 11421900, 11414196, 11421790,11433212, 11433300, 11413100, 11410500)

# Now keep non-duplicate gages (like Bear Sites which have less gaging)
gages_to_drop <- c(11421900, 11414196, 11433212, 11433300, 11413100, 11410500)

combined_df <- combined_df %>% filter(!gageID %in% gages_to_drop)



# GET METADATA ------------------------------------------------------------

# temps: parameterCd = "00010", dailymean = statCd "00003"
# https://nwis.waterdata.usgs.gov/usa/nwis/pmcodes/help?codes_help
# check what daily data is available:
(usgs_daily <- whatNWISdata(siteNumber=combined_df$gageID, service='dv', parameterCd = '00060', statCd='00003') %>% 
   select(site_no, station_nm, dec_lat_va, dec_long_va, huc_cd, 
          data_type_cd, begin_date:count_nu) %>% 
   rename(interval=data_type_cd, huc8=huc_cd, gageID=site_no,
          date_begin=begin_date, date_end=end_date) %>% 
   mutate(yr_begin = year(date_begin),
          yr_end = year(date_end),
          yr_total = yr_end-yr_begin))

# VISUALIZE DATE RANGES FOR SITES -----------------------------------------

# visualize the stations to see date ranges
ggplot() + 
  geom_linerange(data=usgs_daily, 
                 aes(x=forcats::fct_reorder(gageID, huc8), ymin=yr_begin, ymax=yr_end, color=huc8), size=1.1, show.legend = T) + coord_flip() + 
  labs(x="", y="") + scale_color_viridis_d(option = "D")+
  theme_bw(base_family = "Roboto Condensed", base_size = 8) +
ggdark::dark_theme_bw(base_family = "Roboto Condensed", base_size = 8)

#ggsave(filename = "output/figures/usgs_gage_year_ranges_by_interval.pdf",
#       width=11, height = 9.5, units = "in", device = cairo_pdf)


# visualize the stations to see date ranges
ggplot() + 
  geom_linerange(data=usgs_filt, 
                 aes(x=forcats::fct_reorder(gageID, huc8), ymin=yr_begin, ymax=yr_end, color=huc8), size=1.1, show.legend = F) + coord_flip() + 
  labs(x="", y="") +
  scale_color_viridis_d()+
  ggdark::dark_theme_bw(base_family = "Roboto Condensed", base_size = 8)

#ggsave(filename = "output/figures/usgs_gages_data_range_w_8yearsplus_interval.pdf", 
#       width=11, height = 9.5, units = "in", device = cairo_pdf)


# merge datasets
usgs_dataset <- left_join(combined_df, usgs_daily) %>% 
  rename(analysis_yr_begin = year_sta, analysis_yr_end = year_end)

write_csv(usgs_dataset, path = "manuscript/supplemental/Data_S3_usgs_wavelet_analysis.csv")

usgs_dataset %>% select(site:river, analysis_yr_begin:regtype,yr_begin, yr_end) %>% View()

# make spatial
usgs_dataset <- usgs_dataset %>% st_as_sf(coords=c("dec_long_va", "dec_lat_va"), crs=4326)

save(usgs_daily, file="data_output/usgs_wavelet_final_table_sf.rda")

library(mapview)
mapview(usgs_dataset)


# DOWNLOAD DAILY DATA -----------------------------------------------------

# # first filter to only "Daily" stations w > 8yrs data
# usgs_filt <- usgs_daily %>% filter(yr_total > 8)
# usgs_filt_iv <- usgs_event %>% filter(yr_total > 8)
# 
# # Get daily
# usgs_temps_day <- dataRetrieval::readNWISdv(siteNumbers=usgs_filt$site_id, parameterCd = "00010") 
# 
# usgs_temps_day <- usgs_temps_day %>% dataRetrieval::addWaterYear()
# 
# # save out
# save(usgs_temps_day, file = "data/usgs_temps_all_daily_filt8yr.rda")


