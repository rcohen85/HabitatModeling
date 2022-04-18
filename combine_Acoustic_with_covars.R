# This script will combine acoustic presence data at each site and each day with
# covar data from the matching site and day
# One csv file will be created for each species

## SETTINGS --------------------------------------------------------------------
library(tidyverse)
library(stringr)
library(R.matlab)
library(lubridate)

presDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
covarDir = 'J:/Chpt_3/CovarTS'
outDir = 'J:/Chpt_3/ModelData'
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
covarAbbrev = cbind(c("Chl","FSLE","Salinity","SSH","Temperature","VelocityAsp","VelocityMag"),
                    c("Chl","FSLE","Sal","SSH","Temp","VelAsp","VelMag"))
lags = c("Lag7","Lag14","Lag21")

## ACTION ----------------------------------------------------------------------

# Create a list of daily species presence files
presFiles = list.files(presDir, pattern="_Daily.csv", full.names=TRUE)
goodFiles = c(1,2,4,8,11,13:15,17:19)

# find covariate data
varFiles = list.files(covarDir,pattern=".csv",full.names=TRUE)

for (i in goodFiles) {
  
  # Initialize temp data frame
  thisSpecies = numeric()
  
  # Load in a file
  thisFile = data.frame(read.csv(presFiles[i]))
  colnames(thisFile) = c("Date",sites)
  
  # Get species from file name
  species = str_remove(presFiles[i],"_5minBin.csv")
  species = str_remove(species,paste(presDir,'/',sep=""))
  if (str_detect(species,"Atl")){
    species = "Gervais"
  }
  
  # stack presence data from all sites
  thisSpecies = data.frame(Date=as.numeric(rep(as.numeric(as.Date(thisFile[,1],"%d-%b-%Y",origin="1970-01-01")),times=10)),
                           Pres=as.numeric(stack(thisFile[,2:11])[,1]),
                           Site=as.character(rep(sites[1:10],each=dim(thisFile)[1])))

  # Remove NA rows (no-effort days)
  noDat = which(thisSpecies$Pres=="NaN")
  thisSpecies = thisSpecies[-noDat,]
  rownames(thisSpecies) = seq(length=nrow(thisSpecies))
  
  if (i==11) {   # If Risso's, save data for later combination w UD36
    GgTemp = thisSpecies
    next
  } else if (i == 19) {  # If species is UD36, combine w Risso's
    thisSpecies = rbind(GgTemp,thisSpecies)
    thisSpecies = thisSpecies %>% group_by(Date,Site) %>% summarise(Pres=sum(Pres))
    thisSpecies = thisSpecies[order(thisSpecies$Site,thisSpecies$Date),]
  }
  
  for (j in 1:length(varFiles)){ # load covar files one by one
 
    thisVar = data.frame(read.csv(varFiles[j]))
    thisVar[,1] = as.numeric(as.Date(thisVar[,1],origin="1970-01-01"))
    thisVarName = str_remove(str_remove(varFiles[j],paste(covarDir,'/',sep="")),'_TS.csv')
    abbrev = covarAbbrev[which(str_detect(covarAbbrev[,1],thisVarName)),2]
    
    if (abbrev%in%c("Sal","Temp","VelAsp","VelMag")){
      depths = c(0,100,200,300,400,500,600,700)
    } else {
      depths = 0
    }
    
    for (k in 1:length(sites)){ # For each site find data points at this site from all depths and time stamps, plus lags
      
      # match up site & dates
      siteInd = which(!is.na(str_match(thisSpecies$Site,sites[k])))
      putWhere = match(thisVar[,1],thisSpecies$Date[siteInd])
      
      for (l in 1:length(depths)){
        # find covar cols with data for this site & depth (not lag cols!)
        thisSiteDepth = which(!is.na(str_match(colnames(thisVar),
                                               paste(sites[k],as.character(depths[l]),sep=""))) 
                              & !str_detect(colnames(thisVar),"Lag"))
        
        # add data to master data frame
        if (k==1){
        eval(parse(text=paste('thisSpecies$',abbrev,as.character(depths[l]),' = NA',sep="")))}
        eval(parse(text=paste('thisSpecies$',abbrev,as.character(depths[l]),
                              '[siteInd[putWhere[!is.na(putWhere)]]] = thisVar[-which(is.na(putWhere)),thisSiteDepth]',sep="")))
        
        for (m in 1:length(lags)){
          # find covar cols with data for this site & depth & lag
          thisSiteDepth = which(!is.na(str_match(colnames(thisVar),
                                                 paste(sites[k],as.character(depths[l]),sep=""))) 
                                & str_detect(colnames(thisVar),lags[m]))
          
          # add data to master data frame
          if(k==1){
          eval(parse(text=paste('thisSpecies$',abbrev,as.character(depths[l]),lags[m],' = NA',sep="")))}
          eval(parse(text=paste('thisSpecies$',abbrev,as.character(depths[l]),
                                lags[m],'[siteInd[putWhere[!is.na(putWhere)]]] = thisVar[-which(is.na(putWhere)),thisSiteDepth]',sep="")))
        }
      }
    }
  }
  
  thisSpecies = apply(thisSpecies,2,as.character)
  saveName = paste(outDir,'/',species,'_masterDF.csv',sep="")
  write.csv(thisSpecies,saveName,row.names=FALSE)
  
}
