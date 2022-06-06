## SETTINGS -------------------------------------------------------------------
library(lubridate)
library(stringr)
library(ggplot2)
library(viridis)

inDir = 'E:/ModelingCovarData/hycom/0.08deg'
outDir = 'E:/ModelingCovarData/EddyMoviePlots_Temp'
setLon = -74.75
# setLon = c(-74.75,-74.5,-74,-73,-72,-71,-70)
HAT_change = as_date('2017-05-01')
OC_change = as_date('2018-05-01')
HARPs = data.frame(t(data.frame(c(41.06165, -66.35155), # WAT_HZ
                                c(40.26333,-67.9861),  # WAT_OC
                                c(39.83295, -69.98194),  # WAT_NC
                                c(39.19192, -72.22735),  # WAT_BC
                                c(38.37337, -73.36985),  # WAT_WC
                                c(37.16452, -74.46585),  # NFC
                                c(35.30183,-74.8789),  # HAT_A & HAT_B
                                c(33.66992, -75.9977),   # WAT_GS
                                c(32.10527, -77.09067),  # WAT_BP
                                c(30.58295, -77.39002))))  # WAT_BS
colnames(HARPs) = c("lat","lon")
rownames(HARPs) = c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS")

eddy = read.csv("E:/ModelingCovarData/EddyMoviePlots_SSH/AllEddies.csv")
eddy$date = as.Date(eddy$date, format="%Y-%m-%d",tz="UTC")

# fileList = list.files(inDir,pattern = "SSH",full.names = TRUE,recursive=TRUE)
fileList = list.files(inDir,pattern = "Temperature_0",full.names = TRUE,recursive=TRUE)

## CREATE MAPS ----------------------------------------------------------------
# masterData.frontalLat = double()
# masterData.Time = double()

for (i in seq_along(fileList)) {
  # Open a given .Rdata file
  load(fileList[i])
  
  # get 6-digit datestamps from file names
  fileDate = str_extract(fileList[i],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisFileTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  # find if there is a matching date
  matchDate = which(eddy$date==thisFileTime)
  
  if (length(matchDate)==0) {
    next
  } else {
  
  if (thisFileTime>=HAT_change){
    HARPs[7,] = c(35.5841,-74.7499)
  } 
  if (thisFileTime>=OC_change) {
    HARPs[2,] = c(40.22999, -67.97798)
  }
  
  # # Trim the lats
  # trimInd = which(lats>33.5 & lats<38)
  # lats = subset(lats, (lats>33.5 & lats<38))
  # # Trim the SSH data frame
  # data = data[trimInd,]
  
  lons = lons-360

  covar_df = data.frame(data=stack(data.frame(data))[,1],
                        lat=rep(lats,length.out=length(lons)*length(lats)),
                        lon=rep(lons,each=length(lats)))

  # Create eddy data frames
  matchDate = which(eddy$date==thisFileTime)
  thisDateEddy = eddy[matchDate,]
  # Cyclonic data frame
  cyclonic = subset(thisDateEddy,thisDateEddy$polarity=="C")
  cyclonic$centerLon = cyclonic$centerLon-360
  # Anticyclonic data frame
  anticyclonic = subset(thisDateEddy,thisDateEddy$polarity=="A")
  anticyclonic$centerLon = anticyclonic$centerLon-360
  
  # plot for SSH
# map = ggplot(covar_df, aes(x=lon, y=lat)) + geom_tile(aes(fill=data)) +
#     scale_fill_viridis(limits=c(-1.5,1.5)) +
#     geom_point(data=HARPs,aes(x=lon, y=lat), fill="#FDFEFE",color="#FDFEFE",size=2,shape=24) +
#     geom_text(data=cyclonic,aes(x=centerLon, y=centerLat,label="C"), color = "#FF3349", size = 4,fontface="bold")+
#     geom_text(data=anticyclonic,aes(x=centerLon, y=centerLat,label="A"), color = "#EB7A00", size = 4,fontface="bold")+
#     theme(legend.position = "right") +
#     annotate("text",label=as.character(time_temp),x=-79,y=39,color="white",size=6)
  # plot for temp
map = ggplot(covar_df, aes(x=lon, y=lat)) + geom_tile(aes(fill=data)) +
     scale_fill_viridis(limits=c(0,30)) +
     geom_point(data=HARPs,aes(x=lon, y=lat), fill="#FDFEFE",color="#FDFEFE",size=2,shape=24) +
     geom_text(data=cyclonic,aes(x=centerLon, y=centerLat,label="C"), color = "#FF3349", size = 4,fontface="bold")+
     geom_text(data=anticyclonic,aes(x=centerLon, y=centerLat,label="A"), color = "#EB7A00", size = 4,fontface="bold")+
     theme(legend.position = "right") +
     annotate("text",label=as.character(time_temp),x=-79,y=39,color="white",size=6)
  
  # Save this map as a jpeg
  saveName = paste(outDir,'/Eddies_',fileDate,'.jpeg',sep="")
  ggsave(saveName, device="jpeg",width=4,height=2.5,units="in",scale=2.5)
  # mapPlot(coastlineWorld, proj='+proj=moll', col='lightgrey')
  # mapPoints(frontalLon,frontalLat,pch=20,cex=2/3,col='red')
  }
}

# save(,masterData.Time, file=paste(outDir,'frontalLat_TS.Rdata'))


