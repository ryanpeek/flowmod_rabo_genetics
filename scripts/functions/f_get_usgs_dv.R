# function to get data and save/plot:


get_usgs_dv <- function(gage, site, facetYr, filterYrs=F, savedat=T, saveplot=F) {
  gage <- gage
  site <- site
  source("scripts/functions/f_USGS_get_raw_dat.R")
  source("scripts/functions/f_doy.R")
  
  suppressPackageStartupMessages({
    library(tidyverse);
    library(lubridate);
    library(viridis);
    library(ggforce)
  })
  
  # download data:
  get.USGS.raw(gage=gage, river = "dat", sdate = "1900-10-01", 
               saveRaw = F, daily = T)
  
  # filter and format
  dat_dv <- dat_dv %>% filter(!is.na(date)) %>% 
    mutate(flow_cfs=as.numeric(flow_cfs)) %>%
    filter(flow_cfs > 0) %>% 
    add_WYD(., "date")

  print(summary(dat_dv))
  
  # set years for save title 
  maxYr <- max(year(dat_dv$date), na.rm = T)
  minYr <- min(year(dat_dv$date), na.rm= T)

  ## PLOT
  
  facet_yr <- facetYr
  
  vir_D<-ggplot() + 
      geom_line(data=dat_dv, aes(x=date, y=flow_cfs, color=DOWY)) +
      scale_y_log10(labels = scales::comma) + 
      #ylim(0,1000)+
      theme_bw(base_family = "Roboto Condensed")+
      ylab("Discharge (cfs)") + xlab("") +
      ggtitle(paste0(site,": USGS ", 
                     gage, " (", length(unique((dat_dv$WY)))," years)")) +
      scale_color_viridis_c(option = "D") +
      ggforce::facet_zoom(x = WY == facet_yr, shrink=T)

  print(vir_D)
  
  if(saveplot){
    ggsave(filename = paste0("figs/",site,"_",gage,"_",facet_yr,
                             "_usgs_log_flow_facet_zoom.png"), 
           width = 9, height = 7, units = "in", dpi = 300)
  } else{
    print("plot not saved")
  }
  
  if(filterYrs){
    
    summary(dat_dv)
    
    cat("\n","Enter Filter Min Year (YYYY)?","\n\n") # prompt 
       y1<-scan(what="integer",n=1)
    cat("\n","Enter Filter Max Year (YYYY)?","\n\n") # prompt 
       y2<-scan(what="integer",n=1)
    dat_dv <- dat_dv %>% 
      filter(WY > as.integer(y1), WY < as.integer(y2))
    
    summary(dat_dv)
    
    vir_D<-ggplot() + 
      geom_line(data=dat_dv, aes(x=date, y=flow_cfs, color=DOWY)) +
      scale_y_log10(labels = scales::comma) + 
      #ylim(0,1000)+
      theme_bw(base_family = "Roboto Condensed")+
      ylab("Discharge (cfs)") + xlab("") +
      ggtitle(paste0(site,": USGS ", 
                     gage, " (", length(unique((dat_dv$WY)))," years)")) +
      scale_color_viridis_c(option = "D") +
      ggforce::facet_zoom(x = WY == facet_yr, shrink=T)
    
    print(vir_D)
    
    # reset max/min years
    maxYr <- max(year(dat_dv$date), na.rm = T)
    minYr <- min(year(dat_dv$date), na.rm= T)
    
  }
  
  # SAVE OUT: 
  
  # assign to local envir
  assign(x = paste0(site,"_dv"), value=dat_dv)
  
  # create filename
  filename <- paste0(site,"_dv")
  
  # load function to allow naming of save file  
  saveit <- function(..., file) {
    x <- list(...)
    save(list=names(x), file=file, envir=list2env(x))
  }
  
  # if save then save
  if(savedat){
    
    saveit(filename=dat_dv, 
           file = paste0("data/flows/",site,"_dv_USGS_",gage,"_",minYr,"_",maxYr,".rda"))
    
  } else {
    print(paste0(site,"_dv"," in Environment, not saved"))
  }
  
  
}
