# This script will combine the U and V vectors for water velocity
# Water velocity data downloaded from HYCOM, modified for ShinyApp
# Gives magnitude and direction of resultant vector
# LMB 4/4/22

## SETTINGS --------------------------------------------------------------------
library(lubridate)
library(stringr)
library(R.utils)

Indir = 'J:/Chpt_3/HYCOM/0.08deg'
outDir = 'J:/Chpt_3/CovarTS'
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
# depths = c(0,100,200,300,400,500,600,700,800)
depths = c(seq(0,100,by=10),seq(150,300,by=50),seq(400,800,by=100))

# waterU is eastward seawater velocity
# waterV is northward seawater velocity

for (k in 1:length(depths)){
  
  fileListU = list.files(Indir,pattern = paste('U_Velocity_',as.character(depths[k]),'_',sep=""),full.names = TRUE,recursive=FALSE)
  fileListV = list.files(Indir,pattern=paste('V_Velocity_',as.character(depths[k]),'_',sep=""),full.names=TRUE, recursive=FALSE)
  
  # Only change these if using different sites
  OC_change = as_date('2018-05-01') # account for change in
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
  HARPs$Lon1 = HARPs$Lon1+360
  HARPs$Lon2 = HARPs$Lon2+360
  
  
  ## ACTION ----------------------------------------------------------------------
  
  # initialize arrays to hold data from all files matching sites
  masterData.Mag = double()
  masterData.Asp = double()
  masterData.EKE = double()
  masterData.Lat = double()
  masterData.Lon = double()
  masterData.Time = double()
  
  # Load in U-velocity file and V-velocity file from the same date
  
  for (i in seq_along(fileListU)){ # for each file in FileListU
    
    # load the waterU file
    load(fileListU[i])
    # rename waterU files, remove old names from environment
    # this step is necessary to differentiate from water V in calculations
    dataU = data
    rm(data)
    dataU[dataU==0] = NA
    
    # find the matching file at that date in FileListV
    # extract date string from loaded waterU file
    fileDate = str_extract(fileListU[i],"\\d\\d\\d\\d\\d\\d\\d\\d") 
    time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                      str_sub(fileDate,start=5L,end=6L),'-',
                      str_sub(fileDate,start=7L,end=8L),sep="")
    
    thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
    # find the waterV file with the matching date string
    matchDate = which(!is.na(str_match(fileListV, fileDate)))
    
    # load and rename the matching waterV file, if there is one
    if (length(matchDate)>0) {
      
      # load the matching file
      load(fileListV[matchDate])
      # rename waterV files, remove old names from environment
      dataV = data
      rm(data)
      dataV[dataV==0] =NA
      
      # Combine vectors
      # calculate magnitude of resultant vector
      velocityMag = sqrt((dataU*dataU)+(dataV*dataV))
      # calculate angle of resultant vector
      velocityAsp = (pi/2) - atan2(dataV,dataU)
      for (j in 1:dim(velocityAsp)[2]){
        neg = which(velocityAsp[,j]<0)
        missing = which(is.na(velocityAsp[,j]))
        goodIdx = setdiff(neg,missing)
        velocityAsp[goodIdx,j] = velocityAsp[goodIdx,j] + (2*pi)
      }
      # convert from radians to degrees
      velocityAsp = (velocityAsp*180)/pi
      
      # calculate eddy kinetic energy (m^2/s^2)
      # EKE_c = 0.5*(((dataU*100)^2)+((dataV*100)^2))
      EKE = 0.5*((dataU^2)+(dataV^2))
      
      # For Shiny
      # save magnitude, angle, lats, and lons in new file, one per day

      data = velocityMag
      save(data,lats,lons,
           file=paste(Indir,'/','VelocityMag_',as.character(depths[k]),'_',fileDate,'.Rdata',sep=""))
      data = velocityAsp
      save(data,lats,lons,
           file=paste(Indir,'/','VelocityAsp_',as.character(depths[k]),'_',fileDate,'.Rdata',sep=""))
      data = EKE
      save(data,lats,lons,
           file=paste(Indir,'/','EKE_',as.character(depths[k]),'_',fileDate,'.Rdata',sep=""))
      
      
      # For TS
      # initialize arrays to hold relevant data points from this file
      thisFileMag = matrix(nrow=11,ncol=1)
      thisFileAsp = matrix(nrow=11,ncol=1)
      thisFileEKE = matrix(nrow=11,ncol=1)
      thisFileLat = matrix(nrow=11,ncol=1)
      thisFileLon = matrix(nrow=11,ncol=1)
      
      for (m in 1:nrow(HARPs)){ # for each HARP site
        # find data points nearest this HARP site
        if (m==7){ # at HAT, pull points first from site A, then from site B
          if (matchDate<HAT_change){
            sitelat = which.min(abs(HARPs[m,1]-lats))
            sitelon = which.min(abs(HARPs[m,2]-lons))
          } else {
            sitelat = which.min(abs(HARPs[m,3]-lats))
            sitelon = which.min(abs(HARPs[m,4]-lons))
          }
        } else if (m==2) {
          if (matchDate>OC_change){
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
        thisFileMag[m,1] = velocityMag[sitelat-1,sitelon+1]
        thisFileAsp[m,1] = velocityAsp[sitelat-1,sitelon+1]
        thisFileEKE[m,1] = EKE[sitelat-1,sitelon+1]
        thisFileLat[m] = lats[sitelat-1]
        thisFileLon[m] = lons[sitelon+1]
        
      }
      
      # add data points from all files to master data frame
      masterData.Mag = cbind(masterData.Mag, thisFileMag)
      masterData.Asp = cbind(masterData.Asp, thisFileAsp)
      masterData.EKE = cbind(masterData.EKE, thisFileEKE)
      masterData.Lat = cbind(masterData.Lat,thisFileLat)
      masterData.Lon = cbind(masterData.Lon,thisFileLon)
      masterData.Time = cbind(masterData.Time,thisTime)
      
    }
    
    else if (length(matchDate)==0) {
      print(paste('No matching file for waterV found for date: ',as.character(fileDate),sep=""))
      next
    }
    
  }
  
  # Save the time series
  masterData.Data = masterData.Mag
  save(masterData.Data,masterData.Lat,masterData.Lon,masterData.Time,
       file=paste(outDir,'/','VelocityMag','_',as.character(depths[k]),'_TS_ES.Rdata',sep="")) 
  masterData.Data = masterData.Asp
  save(masterData.Data,masterData.Lat,masterData.Lon,masterData.Time,
       file=paste(outDir,'/','VelocityAsp','_',as.character(depths[k]),'_TS_ES.Rdata',sep="")) 
  masterData.Data = masterData.EKE
  save(masterData.Data,masterData.Lat,masterData.Lon,masterData.Time,
       file=paste(outDir,'/','EKE','_',as.character(depths[k]),'_TS_ES.Rdata',sep=""))
}