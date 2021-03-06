# This code will calculate the distance to the gulf stream, saving a time series
# of Gulf Stream frontal position and distance to front for each HARP site

## Settings --------------------------------------------------------------------
library(lubridate)
library(geodist)
library(stringr)

inDir = 'J:/Chpt_3/HYCOM/0.08deg'
outDir = 'J:/Chpt_3/HYCOM/0.08deg'
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
# setLon = -74.75
setLon = -73


#HARP Sites
OC_change = as_date('2018-05-01') # account for change in
HAT_change = as_date('2017-05-01') # account for change in HAT location from site A to B
HARPs = data.frame(t(data.frame(c(41.06165, -66.35155), # WAT_HZ
                                c(40.26333,-67.9861,40.26333,-67.9861, 40.22999, -67.97798),  # WAT_OC
                                c(39.83295, -69.98194),  # WAT_NC
                                c(39.19192, -72.22735),  # WAT_BC
                                c(38.37337, -73.36985),  # WAT_WC
                                c(37.16452, -74.46585),  # NFC
                                c(35.30183,-74.8789,35.5841,-74.7499, 35.5841,-74.7499),  # HAT_A & HAT_B
                                c(33.66992, -75.9977),   # WAT_GS
                                c(32.10527, -77.09067),  # WAT_BP
                                c(30.58295, -77.39002),  # WAT_BS
                                c(30.27818, -80.22085))))  # JAX_D
rownames(HARPs) = sites
colnames(HARPs) = c("Lat1", "Lon1", "Lat2", "Lon2","Lat3","Lon3")

fileList = list.files(inDir,pattern = "SSH",full.names = TRUE,recursive=FALSE)

## Action ----------------------------------------------------------------------

# initialize data frames for saving later
masterData.Frontal = double()
masterData.GeoDist = double()
masterData.Time = double()

# load in SSH files
for (i in seq_along(fileList)){   #for each file in fileList
  
  # Open a given SSH .Rdata file
  load(fileList[i])
  
  # get 6-digit datestamps from file names
  fileDate = str_extract(fileList[i],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisFileTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  # Trim the lats ( 33.5:38 at -74.75, 35:38 at -73)
  trimInd = which(lats>35 & lats<38)
  lats = subset(lats, (lats>35 & lats<38))
  # Trim the SSH data frame
  data = data[trimInd,]
  
  
  # Find the column closest to setLon
  lons = lons-360
  # Adjust lats to remove those that are over land
  # 160 is for -75 setLon
  # lats = lats[1:160]
  # data = data[1:160,]
  colInd = which.min(abs(lons-setLon))
  closeLon = lons[colInd]
  
  # Calculate the first difference of SSH values at column matching closeLon
  minDiffInd = which.min(diff(data[,colInd], lag = 1))
  # Make plots of first difference for each file
  # firstDiff = diff(data[,colInd], lag=1)
  # plot(firstDiff)
  
  # Get the corresponding latitude (gives the last point still in GS)
  minDiffLat = lats[minDiffInd]
  # Add to master data frame
  masterData.Frontal = cbind(masterData.Frontal, minDiffLat)
  
  # Calculate the geodesic distance from HARP site to lat/lon of maxDiff
  if (thisFileTime<HAT_change) {
    frontalLats = rep(minDiffLat, times = 11)
    frontalLons = rep(closeLon, times = 11)
    frontalCoord = t(rbind(frontalLats, frontalLons))
    colnames(frontalCoord) = c("latitude", 'longitude')
    HARPcoords = cbind(HARPs$Lat1, HARPs$Lon1)
    colnames(HARPcoords) = c("latitude", "longitude")
    # calculate geodesic distance in m
    geoDist = geodist(frontalCoord, HARPcoords, paired=TRUE, measure="geodesic")
    geoDist = geoDist/1000
    latOffset = HARPcoords[,1]-frontalCoord[,1]
    southInd = which(latOffset<0)
    geoDist[southInd] = geoDist[southInd]*(-1)
  }
  if (thisFileTime>HAT_change & thisFileTime<OC_change) {
    frontalLats = rep(minDiffLat, times = 11)
    frontalLons = rep(closeLon, times = 11)
    frontalCoord = t(rbind(frontalLats, frontalLons))
    colnames(frontalCoord) = c("latitude", 'longitude')
    HARPcoords = cbind(HARPs$Lat2, HARPs$Lon2)
    colnames(HARPcoords) = c("latitude", "longitude")
    # calculate geodesic distance in m
    geoDist = geodist(frontalCoord, HARPcoords, paired=TRUE, measure="geodesic")
    geoDist = geoDist/1000
    latOffset = HARPcoords[,1]-frontalCoord[,1]
    southInd = which(latOffset<0)
    geoDist[southInd] = geoDist[southInd]*(-1)
  }
  if (thisFileTime>HAT_change & thisFileTime>OC_change) {
    frontalLats = rep(minDiffLat, times = 11)
    frontalLons = rep(closeLon, times = 11)
    frontalCoord = t(rbind(frontalLats, frontalLons))
    colnames(frontalCoord) = c("latitude", 'longitude')
    HARPcoords = cbind(HARPs$Lat3, HARPs$Lon3)
    colnames(HARPcoords) = c("latitude", "longitude")
    # calculate geodesic distance in m
    geoDist = geodist(frontalCoord, HARPcoords, paired=TRUE, measure="geodesic")
    geoDist = geoDist/1000
    latOffset = HARPcoords[,1]-frontalCoord[,1]
    southInd = which(latOffset<0)
    geoDist[southInd] = geoDist[southInd]*(-1)
  }
  # Add geodist to master data frame
  masterData.GeoDist = cbind(masterData.GeoDist, geoDist)
  
  # Add times to master data frame
  masterData.Time = cbind(masterData.Time, thisFileTime)
  
}

# Save TS
save(masterData.GeoDist,masterData.Frontal,masterData.Time,
     file=paste(outDir,'/','GeoDist_',setLon,'_TS.Rdata',sep=""))
# 
# # PLOT STUFF
# load('E:/hycom/0.08deg/TS/GeoDist_-74.5_TS.Rdata')
# hist(masterData.Frontal)
# plot(x=masterData.Time, y=masterData.Frontal)
# 
# plot(x=masterData.Time, y=masterData.GeoDist[7,])
# min(masterData.GeoDist[7,])
# max(masterData.GeoDist[7,])
# 
# 
# library(ggplot2)
# library(rlang)
# load('E:/CovarShinyApp/Covars/Velocity_0_20190525.Rdata')
# data.frame = stack()
# thisData = data.frame(data=stack(data.frame(velocityAngle))[,1],
#                       lat=rep(lats,length.out=length(lons)*length(lats)),
#                       lon=rep(lons-360,each=length(lats)))
# ggplot(thisData, aes(x=lon, y=lat))+geom_tile(aes(fill=data))
# plot(thisData$data)
# hist(thisData$data)
# hist(velocityAngle)

