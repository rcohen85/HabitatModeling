# This script will combine the U and V vectors for water velocity
# Water velocity data downloaded from HYCOM, modified for ShinyApp
# Gives magnitude and direction of resultant vector
# LMB 4/4/22

## SETTINGS --------------------------------------------------------------------
library(lubridate)
library(stringr)
library(R.utils)

dir = 'E:/CovarShinyApp/Covars'

# waterU is eastward seawater velocity
# waterV is northward seawater velocity

fileListU = list.files(dir,pattern = "U_Velocity_0_",full.names = TRUE,recursive=TRUE)
fileListV = list.files(dir,pattern="V_Velocity_0_",full.names=TRUE, recursive=TRUE)

## ACTION ----------------------------------------------------------------------

# Load in U-velocity file and V-velocity file from the same date

for (i in seq_along(fileListU)){ # for each file in FileListU
  
  # load the waterU file
  load(fileListU[i])
  # rename waterU files, remove old names from environment
  # this step is necessary to differentiate from water V in calculations
  dataU = data
  rm(data)
  
  # find the matching file at that date in FileListV
  # extract date string from loaded waterU file
  fileDate = str_extract(fileListU[i],"\\d\\d\\d\\d\\d\\d\\d\\d")
  # find the waterV file with the matching date string
  matchDate = which(!is.na(str_match(fileListV, fileDate)))

  # load and rename the matching waterV file, if there is one
  if (length(matchDate)>0) {
    
    # load the matching file
    load(fileListV[matchDate])
    # rename waterV files, remove old names from environment
    dataV = data
    rm(data)
    
    # Combine vectors
    # calculate magnitude of resultant vector
    velocityMag = sqrt((dataU*dataU)+(dataV*dataV))
    # calculate angle of resultant vector
    velocityAngle = atan2(dataU/dataV)
    velocityAngle[velocityAngle<0] = velocityAngle[velocityAngle<0] + (2*pi)
    # convert from radians to degrees
    velocityAngle = (velocityAngle*180)/pi
    
    # save magnitude, angle, lats, and lons in new file
    save(velocityMag,velocityAngle,lats,lons,
         file=paste(dir,'/','Velocity_0_',fileDate,'.Rdata',sep=""))
  }
  
  else if (length(matchDate)==0) {
    print(paste('No matching file for waterV found. File date',fileDate))
    next
  }
  
  
  
}