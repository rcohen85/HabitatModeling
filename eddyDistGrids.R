# Calculate distance to nearest eddy for each lat/lon point in the region
library(geodist)

inDir = 'J:/Chpt_3/Eddies'

# date range of interest
stDt = as.Date('2016-05-01')
edDt = as.Date('2019-04-30')
dateRange = seq.Date(stDt,edDt,by=1)

lons = seq(278,297,by=0.08)
lats = seq(24,44,by=0.08)

gridDF = data.frame(latitude=rep(lats,length.out=length(lons)*length(lats)),
                    longitude=rep(lons,each=length(lats)))

# load eddy data
eddyData = data.frame(read.csv(paste(inDir,'/AllEddies.csv',sep="")))
eddyData$date = as.Date(eddyData$date,'%d-%b-%Y')

for (j in 1:length(dateRange)){
  
  data = numeric()
  
  # find eddies on this date
  thisDate = which(eddyData$date==dateRange[j])
  
  # calculate min dist to anticyclonic (downwelling) eddy
  Aind = which(eddyData$polarity[thisDate]=='A')
  Acoords = cbind(eddyData$centerLat[thisDate[Aind]],eddyData$centerLon[thisDate[Aind]])
  colnames(Acoords) = c("latitude", 'longitude')
  Adist = geodist(gridDF,Acoords, measure="geodesic")
  Adist = Adist/1000
  data = matrix(apply(Adist,MARGIN=1,min),ncol=length(lons),byrow=FALSE)
  saveName = paste(inDir,'/Grids/AEddyDist_',str_remove_all(dateRange[j],'-'),'.Rdata',sep="")
  save(data,lons,lats,file=saveName)
  
  # calculate min dist to cyclonic (upwelling) eddy
  Cind = which(eddyData$polarity[thisDate]=='C')
  Ccoords = cbind(eddyData$centerLat[thisDate[Cind]],eddyData$centerLon[thisDate[Cind]])
  colnames(Ccoords) = c("latitude", 'longitude')
  Cdist = geodist(gridDF,Ccoords, measure="geodesic")
  Cdist = Cdist/1000
  data = matrix(apply(Cdist,MARGIN=1,min),ncol=length(lons),byrow=FALSE)
  saveName = paste(inDir,'/Grids/CEddyDist_',str_remove_all(dateRange[j],'-'),'.Rdata',sep="")
  save(data,lons,lats,file=saveName)

  
}