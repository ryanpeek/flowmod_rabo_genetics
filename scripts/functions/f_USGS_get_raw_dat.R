### Download USGS Gage Data 
### by R. Peek 2017

# All data from USGS NWIS website, in RDB format. 
# Function will download daily or instantaneous (IV) 15-min data

 # For 15 min data: 
  #  http://waterservices.usgs.gov/nwis/iv/?format=rdb&sites=11427000&startDT=2011-10-01&endDT=2012-08-23&parameterCd=00060,00065,00010 ##

 # For for multiple gages in one call "&sites=01646500,06306300"
	#  http://waterservices.usgs.gov/nwis/iv/?format=rdb&sites=11413000,11427000&startDT=2009-10-01&endDT=2015-07-26&parameterCd=00060,00065,00010

# INSTANTANEOUS -----------------------------------------------------------

get.USGS.raw<-function(gage, river, sdate="2010-10-01", edate=Sys.Date(),
                   saveRaw=FALSE, daily=FALSE){
  ##  SITES
  ##  11413000 NF Yuba,  Goodyears Bar
  ##  11427000 NF American,  Clementine Dam
  ##  11285500 Tuo R at Wards Ferry
  
  ## These may not work...
  ##  11210100 SF Kaweah,   Three Rivers
  ##  11264500 Merced, 	Happy Isles Bridge
  ##  11276600 Tuolumne, Early Intake
  ##  11335000 Cosumnes,  Michigan Bar  
  ##  11401500 Feather,  Indian Creek nr Crescent Mill

  # load packages
  if(!require(lubridate)) { install.packages("lubridate"); library(lubridate)}
  if(!require(data.table)) { install.packages("data.table"); library(data.table)}
  if(!require(dplyr)) { install.packages("dplyr"); library(dplyr, warn.conflicts = F)}
  
  # set basic variables
  gage = gage; sdate = sdate; river = river; edate=edate
  
  if(daily){
    timeperiod <- "dv"
  } else {
    timeperiod <- "iv"
  }
  
  # read raw data
  rawdat<-fread(input = paste("http://waterservices.usgs.gov/nwis/",
                              timeperiod,"/?format=rdb&sites=",
                              gage,"&startDT=",sdate,"&endDT=",edate,
                              "&parameterCd=00060,00065",sep=""), 
                sep="\t", skip = 28) 
  
  # test line: paste("http://waterservices.usgs.gov/nwis/iv/?format=rdb&sites=",11427000,"&startDT=","2010-10-01","&endDT=","2015-07-30","&parameterCd=00060,00065",sep="")
  
  riv.df <- as.data.frame(rawdat) # read data
  
  # adjust/format the datetime to POSIXct (good for plotting w ggplot)
  # to check time zones, use:   grep("^America", OlsonNames(),value=TRUE)
  
  if(daily){
    riv.df[[3]]<-ymd(riv.df[[3]])
    riv <-data.frame(riv.df) # subset to columns of interest
    colnames(riv)<- c("agency","gageNo", "date", "flow_cfs","cd_00060") # rename cols
    riv[,5]<-as.factor(riv[,5]) # this is not totally necessary
  } else {
    
    riv.df$datetime<-ymd_hm(paste(riv.df[[3]]," ",riv.df[[4]],sep=""),tz="America/Los_Angeles")
    riv <-data.frame(riv.df) 
    colnames(riv)<- c("agency","gageNo", "datetime_raw","TZ", "flow_cfs","cd_00060","stage_ft", "cd_00065", "datetime") # rename cols
    riv[,6]<-as.factor(riv[,6])
    riv[,8]<-as.factor(riv[,8])
    riv <- select(riv, 1:2,9,5:8) # simplify cols (reg rid of datetime_raw and TZ)
  }
  
  print(summary(riv))
  print(head(riv))
  
  # Save data ---------------------------------------------------------------   
  # assumes data can be saved in a "data" folder in a root dir
  
  if(!require(readr)) { install.packages("readr"); library(readr)}

  if(saveRaw){
    write_rds(riv, path=paste0("data/",river,"_",sdate,"_",edate, "_",timeperiod,"_USGS.rds"),compress = "gz")
    print(paste("data saved here: ",getwd(),sep=""))
  } else {
    cat("USGS raw data not saved.\n")
  }
  
  # save data to workspace/environment
  
  river1<-paste0(river,"_",timeperiod)
  assign(river1,riv, envir=.GlobalEnv) # print to workspace
  
  cat("All finished... data in current environment \n")
}  
