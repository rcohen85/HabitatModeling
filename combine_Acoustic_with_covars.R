# This script will combine acoustic presence data at each site and each day with
# covar data from the matching site and day
# One csv file will be created for each species

## SETTINGS --------------------------------------------------------------------
library(stringr)
library(R.matlab)

inDir = 'E:/NEWdailyTotals'
outDir = 'E:/ModelingCovarData'
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')

## ACTION ----------------------------------------------------------------------

# Create a list of daily totals files
dailyTotsFiles = list.files(inDir, pattern=".mat", full.names=TRUE)

# Initialize data frames for each species
master.Blainville = double()
master.UD26 = double()
master.UD28 = double()
master.Cuvier = double()
master.Gervais = double()
master.Kogia = double()
master.Risso = double()
master.Sowerby = double()
master.SpermWhale = double()
master.True = double()

masterDfList = c("master.Blainville","master.UD26","master.UD28","master.Cuvier","master.Gervais","master.Kogia","master.Risso","master.Sowerby","master.SpermWhale","master.True")

for (i in seq_along(dailyTotsFiles)) {
  
  # Initialize temp data frames for each species
  tempBlainville = double()
  tempUD26 = double()
  tempUD28 = double()
  tempCuvier = double()
  tempGervais = double()
  tempKogia = double()
  tempRisso_UD36 = double()
  tempSowerby = double()
  tempSpermWhale = double()
  tempTrue = double()
  
  dfList = c("tempBlainville","tempUD26","tempUD28","tempCuvier","tempGervais","tempKogia","tempRisso_UD36","tempSowerby","tempSpermWhale","tempTrue")
  speciesInd = data.frame(c(dfList),c(1,4,5,9,21,13,3,18,19,20),c(1,4,5,9,21,13,16,18,19,20))
  colnames(speciesInd) = c("Species","Ind1","Ind2")
  
  # Load in a file
  thisDailyTotsFile = readMat(dailyTotsFiles[i])
  
  # Get species list from file
  speciesList = t(data.frame(thisDailyTotsFile$spNameList))
  
  # Get site from file name
  for (j in seq_along(sites)) {
    if (isTRUE(str_detect(dailyTotsFiles[i], sites[j]))) {
      site = sites[j]
    }
  }
  
  # Get daily totals from file
  dailyTots = thisDailyTotsFile$dailyTots
  
  # Get date from daily totals
  matFileDate = dailyTots[,1]
  # Convert to R dates
  fileDate = as.POSIXct((matFileDate-719529)*86400,format='%Y-%m-%d',origin='1970-01-01',tz="UTC")
  # fileDate = as.numeric(fileDate)
  
  # Get species list from file
  speciesList = t(data.frame(thisDailyTotsFile$spNameList))
  
  
  
  # Add site and dates to data frames
  for (j in seq_along(dfList)) {
    # Add fileDate as a column in each data frame
    eval(parse(text=paste(dfList[j],"=cbind(",dfList[j],",as.Date(fileDate))")))
    # Add site as a column in each data frame
    eval(parse(text=paste(dfList[j],"=cbind(",dfList[j],",site)")))
    # Make them all data frames
    eval(parse(text=paste(dfList[j],"=data.frame(",dfList[j],")")))
    # Change column names
    eval(parse(text=paste("colnames(",dfList[j],")=c('Date','Site')")))
   
    if (j!=7) {   # If species is not Risso's
      # Find the daily totals column for this species
      thisSpeciesDT = dailyTots[,speciesInd$Ind1[j]+1]
      # Add it to the data frames
      eval(parse(text=paste(dfList[j],"=cbind(",dfList[j],",thisSpeciesDT)")))
    }
    if (j == 7) {
      thisSpeciesDT1 = dailyTots[,speciesInd$Ind1[j]+1]
      thisSpeciesDT2 = dailyTots[,speciesInd$Ind2[j]+1]
      thisSpeciesDT = c(thisSpeciesDT1,thisSpeciesDT2)
      tempRisso_UD36 = data.frame(Date=rep(tempRisso_UD36$Date, length.out=(length(tempRisso_UD36$Date))*2),
                              Site=rep(tempRisso_UD36$Site, length.out=(length(tempRisso_UD36$Site))*2),
                              thisSpeciesDT=thisSpeciesDT)
    }
  }
  
 for (j in seq_along(masterDfList)) {
   eval(parse(text=paste(masterDfList[j],"=rbind(",masterDfList[j],",",dfList[j],")")))
 }
  
}
