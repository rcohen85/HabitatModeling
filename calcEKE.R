# Calculate EKE from u- and v-velocities
# EKE calculation should be based on u' and v', which are the surface velocity fluctuations,
# obtained by subtracting the annual means (u_bar and v_bar) from the total u and v
library(zoo)
library(stringr)

inDir = 'J:/Chpt_3/HYCOM/0.08deg'
outDir = 'J:/Chpt_3/CovarTS'
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
depths = c(seq(0,100,by=10),seq(150,300,by=50),seq(400,800,by=100))

# Only change these if using different sites
OC_change = as.Date('2018-05-01') # account for change in
HAT_change = as.Date('2017-05-01') # account for change in HAT location from site A to B
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
HARPs$Lon1 = HARPs$Lon1+360
HARPs$Lon2 = HARPs$Lon2+360

for (k in 1:length(depths)){

  fileListU = list.files(Indir,pattern = paste('U_Velocity_',as.character(depths[k]),'_',sep=""),full.names = TRUE,recursive=FALSE)
  fileListV = list.files(Indir,pattern=paste('V_Velocity_',as.character(depths[k]),'_',sep=""),full.names=TRUE, recursive=FALSE)

  masterData.Data = double()
  masterData.Lat = double()
  masterData.Lon = double()
  masterData.Time = double()
  allU = double()
  allV = double()
  depthDate = double()

  for (i in seq_along(fileListU)){ # for each file in FileListU

    # extract date string from loaded waterU file
    fileDate = str_extract(fileListU[i],"\\d\\d\\d\\d\\d\\d\\d\\d")
    time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                      str_sub(fileDate,start=5L,end=6L),'-',
                      str_sub(fileDate,start=7L,end=8L),sep="")

    thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")

    # find the waterV file with the matching date string
    matchDate = which(!is.na(str_match(fileListV, fileDate)))

    if (length(matchDate)>0) { # if there are both U and V files for this date, load them

      # load the waterU file
      load(fileListU[i])
      dataU = stack(data.frame(data))[,1]
      dataU[dataU==0] = NA
      allU = rbind(allU,dataU)

      load(fileListV[matchDate])
      dataV = stack(data.frame(data))[,1]
      dataV[dataV==0] =NA
      allV = rbind(allV,dataV)

      depthDate = c(depthDate,thisTime)
    }
  }

  # Calculate moving average of velocities (to remove very short-term fluctuations)
  smoothU = apply(allU,MARGIN=2,rollmean,k=5,fill=NA)
  smoothV = apply(allV,MARGIN=2,rollmean,k=5,fill=NA)
  
  # Calculate yearly means and subtract from averages to get u' and v'
  ints = seq.Date(as.Date("2016-05-01"),as.Date("2019-05-01"),by=365)
  uPrime = numeric()
  vPrime = numeric()
  gridDates = numeric()
  for (i in 1:(length(ints)-1)){
    
    # find grids in this year
    whichRows = which(depthDate>=ints[i] & depthDate<ints[i+1])
    gridDates = c(gridDates,depthDate[whichRows])
    
    # average each grid point across the year 
    uMean = apply(smoothU[whichRows,],MARGIN=2,mean,na.rm=TRUE)
    vMean = apply(smoothV[whichRows,],MARGIN=2,mean,na.rm=TRUE)
    
    # subtract yearly mean from each day
    uFluct = t(t(smoothU[whichRows,]) - uMean) 
    vFluct = t(t(smoothV[whichRows,]) - vMean)
    
    uPrime = rbind(uPrime,uFluct)
    vPrime = rbind(vPrime,vFluct)
  }
  
  # Calculate EKE from u' and v'
  EKE = 0.5*(((vPrime*100)^2) + ((vPrime*100)^2))
  
  # Save daily grids and create time series at each HARP site
  for (i in 1:dim(EKE)[1]){
    
    printDate = str_remove_all(as.Date(gridDates[i]),'-')
    
    data = matrix(EKE[i,],ncol=length(lons),byrow=FALSE)
    saveName = paste(inDir,'/','EKE_',as.character(depths[k]),'_',printDate,'.Rdata',sep="")
    save(data,lats,lons,file=saveName)
    
    thisDateEKE = matrix(nrow=11,ncol=1)
    for (m in 1:nrow(HARPs)){ # for each HARP site
      # find data points nearest this HARP site
      if (m==7){ # at HAT, pull points first from site A, then from site B
        if (gridDates[i]<HAT_change){
          sitelat = which.min(abs(HARPs[m,1]-lats))
          sitelon = which.min(abs(HARPs[m,2]-lons))
        } else {
          sitelat = which.min(abs(HARPs[m,3]-lats))
          sitelon = which.min(abs(HARPs[m,4]-lons))
        }
      } else if (m==2) {
        if (gridDates[i]>OC_change){
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
      thisDateEKE[m,1] = mean(data[(sitelat-2):(sitelat-1),(sitelon+1):(sitelon+2)],na.rm=TRUE)
    }
    
    masterData.Data = cbind(masterData.Data, thisDateEKE)
    masterData.Lat = cbind(masterData.Lat,lats)
    masterData.Lon = cbind(masterData.Lon,lons)
    masterData.Time = cbind(masterData.Time,gridDates[i])
    
  }
  
  save(masterData.Data,masterData.Lat,masterData.Lon,masterData.Time,
       file=paste(outDir,'/','EKE','_',as.character(depths[k]),'_TS_ES_new.Rdata',sep=""))
}