# Calculate distance to nearest cyclonic/anticyclonic eddy for a set of HARP sites
library(geodist)

inDir = 'J:/Chpt_3/Eddies'
outDir = 'J:/Chpt_3/CovarTS'

# date range of interest
stDt = as.Date('2016-05-01')
edDt = as.Date('2019-04-30')
dateRange = seq.Date(stDt,edDt,by=1)

# HARP sites
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
OC_change = as.Date('2018-05-01') # account for change in
HAT_change = as.Date('2017-05-01') # account for change in HAT location from site A to B
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


## --------------------------------------------

# load eddy data
eddyData = data.frame(read.csv(paste(inDir,'/AllEddies.csv',sep="")))
eddyData$date = as.Date(eddyData$date,'%d-%b-%Y')
eddyData$centerLon = eddyData$centerLon-360

masterData.eddyDistA = double()
masterData.eddyDistC = double()
masterData.Time = double()

for (j in dateRange){
  
  # save time stamp
  masterData.Time = cbind(masterData.Time, j)
  
  # find eddies on this date
  thisDate = which(eddyData$date==j)
  
  # get correct HARP coordinates based on the date
  if (j<HAT_change){
    HARPcoords = HARPs[,1:2]
  } else if (j>=HAT_change & j<OC_change){
    HARPcoords = HARPs[,3:4]
  } else if (j>=OC_change) {
    HARPcoords = HARPs[,5:6]}
  
  colnames(HARPcoords) = c("latitude", 'longitude')
  
  # calculate min dist to anticyclonic (downwelling) eddy
  Aind = which(eddyData$polarity[thisDate]=='A')
  Acoords = cbind(eddyData$centerLat[thisDate[Aind]],eddyData$centerLon[thisDate[Aind]])
  colnames(Acoords) = c("latitude", 'longitude')
  Adist = geodist(HARPcoords,Acoords, measure="geodesic")
  Adist = Adist/1000
  masterData.eddyDistA = cbind(masterData.eddyDistA,apply(Adist,MARGIN=1,min))
  
  # calculate min dist to cyclonic (upwelling) eddy
  Cind = which(eddyData$polarity[thisDate]=='C')
  Ccoords = cbind(eddyData$centerLat[thisDate[Cind]],eddyData$centerLon[thisDate[Cind]])
  colnames(Ccoords) = c("latitude", 'longitude')
  Cdist = geodist(HARPcoords,Ccoords, measure="geodesic")
  Cdist = Cdist/1000
  masterData.eddyDistC = cbind(masterData.eddyDistC,apply(Cdist,MARGIN=1,min))
  
}

# Save TS
save(masterData.eddyDistC,masterData.eddyDistA,masterData.Time,
     file=paste(outDir,'/','EddyDist_TS.Rdata',sep=""))




  