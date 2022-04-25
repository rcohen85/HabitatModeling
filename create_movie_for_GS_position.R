# This script will create a map of our calculated GS position for each day
# The maps are saved into a new directory
# The maps are compiled into a movie
# LMB 4/4/22

## SETTINGS -------------------------------------------------------------------
library(lubridate)
library(stringr)
library(ggplot2)
library(viridis)

inDir = 'E:/hycom/0.08deg'
outDir = 'E:/hycom/0.08deg/movie/'
setLon = -74.75
# setLon = c(-74.75,-74.5,-74,-73,-72,-71,-70)
HAT_change = as_date('2017-05-01')
HARPs = data.frame(t(data.frame(c(38.37337, -73.36985),  # WAT_WC
                                c(37.16452, -74.46585),  # NFC
                                c(35.30183,-74.8789),  # HAT_A & HAT_B
                                c(33.66992, -75.9977))))  # WAT_GS
colnames(HARPs) = c("lat","lon")
rownames(HARPs) = c("WC","NFC","HAT","GS")

fileList = list.files(inDir,pattern = "SSH",full.names = TRUE,recursive=TRUE)

## CREATE MAPS ----------------------------------------------------------------
# masterData.frontalLat = double()
# masterData.Time = double()

for (i in seq_along(fileList)) {
  # Open a given SSH .Rdata file
  load(fileList[i])
  
  # get 6-digit datestamps from file names
  fileDate = str_extract(fileList[i],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisFileTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  if (thisFileTime>=HAT_change){
    HARPs[3,] = c(35.5841,-74.7499)
  } 
  
  # Trim the lats
  trimInd = which(lats>33.5 & lats<38)
  lats = subset(lats, (lats>33.5 & lats<38))
  # Trim the SSH data frame
  data = data[trimInd,]
  
  # Find the column closest to setLon
  lons = lons-360
  # Adjust lats to remove those that are over land
  # 160 is for -75 setLon
  # lats = lats[1:160]
  # data = data[1:160,]
  colInd = as.numeric()
  closeLon = as.numeric()
  minDiffInd = as.numeric()
  minDiffLat = as.numeric()
  frontalLat = as.numeric()
  frontalLon = as.numeric()
  
  for (j in 1:length(setLon)){
    
  colInd = c(colInd,which.min(abs(lons-setLon[j])))
  closeLon = c(closeLon,lons[colInd[j]])
  
  # Calculate the first difference of SSH values at column matching closeLon
  minDiffInd = c(minDiffInd,which.min(diff(data[,colInd[j]], lag = 1)))
  
  # # Calculate max fsle value at column matching closeLon
  # minDiffInd = c(minDiffInd,which.min(data[,colInd[j]]))
  
  # Calculate the frontal position of GS (last point still in the GS)
  frontalLat = c(frontalLat,lats[minDiffInd[j]])
  frontalLon = c(frontalLon,closeLon[j])
  
  }
  
  frontal_df = data.frame(lat=frontalLat,lon=frontalLon)
  
  # masterData.frontalLat = cbind(masterData.frontalLat,frontalLat)
  # masterData.Time = cbind(masterData.Time, thisFileTime)

  # Make plots of first difference for each file
  # firstDiff = diff(data[,colInd], lag=1)
  # plot(firstDiff)
  
  # create a data frame for ggplot
  covar_df = data.frame(data=stack(data.frame(data))[,1],
             lat=rep(lats,length.out=length(lons)*length(lats)),
             lon=rep(lons,each=length(lats)))
  # # replace zero and low values with NA
  # covar_df$data[covar_df$data == 0] = NA
  # covar_df$data[covar_df$data<25] = NA
  
  # Plot the frontal position of GS
  # map = ggplot(ssh_df, aes(x=lon, y=lat)) + geom_tile(aes(fill=data)) + 
  #   geom_point(x=frontalLon, y=frontalLat, aes(color = "red", size = 6)) +
  #   scale_fill_viridis()
    
  map = ggplot(covar_df, aes(x=lon, y=lat)) + geom_tile(aes(fill=data)) +
      scale_fill_viridis(limits=c(-1.5,1.5)) +
      geom_point(data=HARPs,aes(x=lon, y=lat), fill="#FDFEFE",color="#FDFEFE",size=2,shape=24) +
      geom_point(data=frontal_df,aes(x=lon, y=lat), color = "#FF3349", size = 2)+
      annotate("text",label=as.character(time_temp),x=-79,y=39,color="white",size=6)

  # Save this map as a jpeg
  saveName = paste(outDir,'frontalPosition_',fileDate,'.jpeg',sep="")
  ggsave(saveName, device="jpeg",width=4,height=2.5,units="in",scale=2.5)
  # mapPlot(coastlineWorld, proj='+proj=moll', col='lightgrey')
  # mapPoints(frontalLon,frontalLat,pch=20,cex=2/3,col='red')
}

# save(,masterData.Time, file=paste(outDir,'frontalLat_TS.Rdata'))

## CREATE ANIMATION ----------------------------------------------------------
# library(imagemagick)
# 
# plotList = list.files(inDir,pattern = ".jpeg",full.names = TRUE,recursive=TRUE)
# 
# animation = image_animate(plotList, fps=2, optimize=2)
