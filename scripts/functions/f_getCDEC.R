#### FUNCTION TO PULL DATA FROM CDEC STATIONS
#### RYAN PEEK, 2016, 
#### CENTER FOR WATERSHED SCIENCES UC DAVIS

get.CDEC<-function( # Function to pull CDEC data
  station,  # Station is 3 letter abbreviation
  duration, # Duration is E=event, D=Daily, H=Hourly
  sensor,   # sensor is number, see below
  start,    # "YYYY-MM-DD"
  end,      # "YYYY-MM-DD"
  csv=F){     # export to CSV?
  
  # List of Real-Time Stations: http://cdec.water.ca.gov/misc/realStations.html
  # List of Daily Stations: http://cdec.water.ca.gov/misc/dailyStations.html 
  # List of sensors:  http://cdec.water.ca.gov/misc/senslist.html
  
  ## STATIONS USED FREQUENTLY
  # SGP: Sugar Pine (AMR)
  # HYS: Huysink (AMR)
  # TGC: Tuolumne Grand Canyon (TUO)
  # DNP: Don Pedro (TUO)
  
  ## SENSORS MOST COMMONLY USED
  # 1  stage (ft)
  # 20 flow cfs
  
  # 2  rain accum (in)
  # 16 precip tippingbucket (in)
  # 45 ppt incremental
  
  # 3  snow water content (in)
  # 18 snow depth (in)
  
  # 6  reservoir elevation (ft)
  # 15 reservoir storage (ac-ft)
  # 76 reservoir inflow
     
  # 25 water temp 
  # 4  air temp 
  # 12 Rel humidity
  # 17 Baro pressure (atmospheric pressure)
  
  ## EXAMPLE URL:  
  # http://cdec.water.ca.gov/cgi-progs/queryCSV?station_id=OXB&dur_code=E&sensor_num=20&start_date=2011/10/01&end_date=2012/09/30
  
  data <- read.table(paste("http://cdec.water.ca.gov/cgi-progs/queryCSV?station_id=", station,"&dur_code=",duration,"&sensor_num=", sensor,"&start_date=", start, "&end_date=",end, "&data_wish=View+CSV+Data", sep=""),
                     header=T, sep=",",quote="'", skip=1, na.strings="m",colClasses = c("character", "character", "numeric"))
  names(data)[1] <- "date"
  names(data)[3]<- "data"
  data$date<-as.Date(data$date,"%Y%m%d")
  
  ## format date and time
  names(data) <- c("date","time", paste("sensor_",sensor,sep=""))
  data$date<-strptime(data$date,format="%Y-%m-%d") # convert to datetime
  data$time<-as.character(data$time)
  data$datetime<-paste(data$date," ",data$time,sep="") # create a datetime column by pasting them together
  data$datetime<-as.POSIXct(strptime(data$datetime,format="%Y-%m-%d %H%M")) # convert to datetime from whatever format
  data$station<-station
  summary(data)
  str(data)
  
  #   ## ask if user wants to change the directory for save purposes    
  #   cat("\n","Use current directory? Y or N","\n\n") # prompt 
  #   y<-scan(what="character",n=1)
  #   ifelse(y=="N",setwd(choose.dir()),getwd())
  
  ## ask if user wants to save to csv or use in dataframe
  #cat("\n","Write file to csv? Y or N","\n\n") # prompt 
  #z<-scan(what="character",n=1)
  if(csv){write.csv(data, file=paste("cdec_", station,"_sensor-",sensor,"_",start,"_to_",end,".csv",sep=""),row.names=FALSE)
             print(paste("file downloaded and saved here: ",getwd(),sep="")) # show message
  } else{
    cat("Output to dataframe only\n")}
  cdec<-paste0(station)
  assign(cdec, data,envir = .GlobalEnv) # print to workspace
  cat("All Finished! Available in current dataframe...\n")
}

