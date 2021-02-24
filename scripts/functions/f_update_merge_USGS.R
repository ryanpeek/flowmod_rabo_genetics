# Update and merge USGS gage data
# R. Peek 2017

# assumes master file is in the "data_output" folder of root dir
library(tidyverse)

update_merge_USGS <- function(
  site, # the name of the site
  #gageID, # the gage ID
  updated_df, # the df from get.USGS.raw
  dat_interval="iv"){ # either "iv" for instantaneous or "dv" for daily 
  
  # get file name for given variables
  dat_File <- list.files(path = "data_output", pattern = paste0(site, "_", dat_interval, "_USGS_updated"))
  
  # get the last updated date: 
  #updatedDate <- gsub("[A-z \\.\\(\\)]", "", dat_File)
  
  # Load Existing Data (15 min, since 1980, Daily since 1930/40)
  dat_archived <- read_rds(paste0("data_output/", dat_File))
  
  # actual rbind and merge (use distinct to remove dups)
  if(dat_interval=="iv"){
    dat_update <- bind_rows(dat_archived, updated_df) %>% 
    distinct(datetime, .keep_all = T)
  } else {
    dat_update <- bind_rows(dat_archived, updated_df) %>% 
      distinct(date, .keep_all = T)
  }
  
  # write data
  write_rds(dat_update, path = paste0("data_output/",site, 
                                      "_", dat_interval, "_USGS_updated_", 
                                      Sys.Date(),".rds"), compress = "gz")
  
  # add name for site
  site_updated_df <- paste0(site, "_", dat_interval, "_updated")
  
  # assign to environment
  assign(site_updated_df, dat_update, envir=.GlobalEnv) # rename with site
  
  # rm temp dataframes
  rm(list = ls(pattern = "^dat_")) # rm old df
  rm(site, site_updated_df)
}


