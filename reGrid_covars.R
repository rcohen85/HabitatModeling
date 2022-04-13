library(EFDR)
library(stringr)
library(lubridate)
library(ncdf4)
library(ggplot2)

inDir = 'J:/Chpt_3/FSLE/0.04deg'
outDir = 'J:/Chpt_3/FSLE/0.08deg'
covar = 'FSLE'
newx = seq(278,297,by=0.08)
newy = seq(24,44,by=0.08)
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
OC_change = as_date('2018-05-01') # account for change in OC site location
HAT_change = as_date('2017-05-01') # account for change in HAT location from site A to B
HARPs = data.frame(t(data.frame(c(41.06165, -66.35155), # WAT_HZ
                                c(40.26333,-67.9861,40.22999, -67.97798),  # WAT_OC
                                c(39.83295, -69.98194),  # WAT_NC
                                c(39.19192, -72.22735),  # WAT_BC
                                c(38.37337, -73.36985),  # WAT_WC
                                c(37.16452, -74.46585),  # NFC
                                c(35.30183,-74.8789,35.5841,-74.7499),  # HAT_A & HAT_B
                                c(33.66992, -75.9977),   # WAT_GS
                                c(32.10527, -77.09067),  # WAT_BP
                                c(30.58295, -77.39002),  # WAT_BS
                                c(30.27818, -80.22085))))  # JAX_D
rownames(HARPs) = sites
colnames(HARPs) = c("Lat1", "Lon1", "Lat2", "Lon2")

## -------------------

fileList = dir(inDir,covar)
HARPs$Lon1 = HARPs$Lon1+360
HARPs$Lon2 = HARPs$Lon2+360

# initialize arrays to hold data from all files matching sites
masterData.Data=double()
masterData.Lat=double()
masterData.Lon=double()
masterData.Time=double()

for (i in 1:length(fileList)){

  load(paste(inDir,'/',fileList[i],sep=""))
  if (typeof(data)=="list"){
    data=data.matrix(data)
  }
  data[is.na(data)] = mean(data,na.rm=TRUE)
  
  # get 6-digit datestamp from file name
  fileDate = str_extract(fileList[i],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  thisData = data.frame(z=stack(data.frame(data))[,1],
                        y=rep(lats,length.out=length(lons)*length(lats)),
                        x=rep(lons,each=length(lats)))
  # ggplot(thisData,aes(x=x,y=y))+geom_tile(aes(fill=z))
  
  newGrid = regrid(thisData,n1=length(newx),n2=length(newy),method="idw")
  # ggplot(newGrid,aes(x=x,y=y))+geom_tile(aes(fill=z))
  data = matrix(newGrid$z,ncol=length(newx),byrow=TRUE)
  lons = sort(unique(newGrid$x))
  lats = sort(unique(newGrid$y))

  save(data,lats,lons,file=paste(outDir,'/',fileList[i],sep=""))
  
  ## Get data at each HARP site
  thisFileData = matrix(nrow=11,ncol=1)
  thisFileLat = matrix(nrow=11,ncol=1)
  thisFileLon = matrix(nrow=11,ncol=1)
  
  for (m in 1:nrow(HARPs)){ # for each HARP site
    
    # find data points nearest this HARP site
    if (m==7){ # at HAT, pull points first from site A, then from site B
      if (thisTime<HAT_change){
        sitelat = which.min(abs(HARPs[m,1]-lats))
        sitelon = which.min(abs(HARPs[m,2]-lons))
      } else {
        sitelat = which.min(abs(HARPs[m,3]-lats))
        sitelon = which.min(abs(HARPs[m,4]-lons))
      }
    } else if (m==2) {# at OC, account for slight change in site location
      if (thisTime<OC_change){
        sitelat = which.min(abs(HARPs[m,1]-lats))
        sitelon = which.min(abs(HARPs[m,2]-lons))
      } else {
        sitelat = which.min(abs(HARPs[m,3]-lats))
        sitelon = which.min(abs(HARPs[m,4]-lons))
      }
    } else {
      sitelat = which.min(abs(HARPs[m,1]-lats))
      sitelon = which.min(abs(HARPs[m,2]-lons))
    }
    
    # grab data values at this HARP site
    thisFileData[m,1] = data[sitelat,sitelon]
    thisFileLat[m] = lats[sitelat]
    thisFileLon[m] = lons[sitelon]
    
  }
  
  masterData.Data = cbind(masterData.Data, thisFileData)
  masterData.Lat = cbind(masterData.Lat,thisFileLat)
  masterData.Lon = cbind(masterData.Lon,thisFileLon)
  masterData.Time = cbind(masterData.Time,thisTime)
  
}

save(masterData.Data,masterData.Lat,masterData.Lon,masterData.Time,
     file=paste(outDir,'/',covar,'_TS.Rdata',sep=""))
