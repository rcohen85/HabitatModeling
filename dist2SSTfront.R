library(raster)


inDir = 'J:/Chpt_3/HYCOM/0.08deg'
fileList = dir(inDir,"Temperature_0")

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

masterData.time=numeric()
masterData.SSTdist=numeric()

for (i in 1:length(fileList)){
  
  load(paste(inDir,"/",fileList[i],sep=""))
  
  # get 6-digit datestamps from file names
  fileDate = str_extract(fileList[i],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  # convert land grid points (zeros) to NAs
  land = which(data==0)
  data[land] = NA
  
  # find edges
  canEd = image_canny_edge_detector((data))
  
  # get rid of land edges
  canEd$edges[land] = NA
  
  # get correct HARP coordinates based on the date
  if (thisTime<HAT_change){
    HARPcoords = HARPs[,1:2]
  } else if (thisTime>=HAT_change & j<OC_change){
    HARPcoords = HARPs[,3:4]
  } else if (thisTime>=OC_change) {
    HARPcoords = HARPs[,5:6]}
  
  colnames(HARPcoords) = c("latitude", 'longitude')
  
  # calculate distance from each HARP site to nearest edge
  
  
}
