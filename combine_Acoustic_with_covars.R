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
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS')
covarAbbrev = cbind(c("Chl","FSLE","Salinity","SSH","Temperature","VelocityAsp","VelocityMag","EKE","AEddyDist","CEddyDist"),
                    c("Chl","FSLE","Sal","SSH","Temp","VelAsp","VelMag","EKE","AEddyDist","CEddyDist"))
# lags = c(7,14,21,28,42,56)
OC_change = as_date('2018-05-01') # account for change in OC site location
HAT_change = as_date('2017-05-01') # account for change in HAT location from site A to B
## ACTION ----------------------------------------------------------------------

# Create a list of daily species presence files
presFiles = list.files(presDir, pattern="_Daily.csv", full.names=TRUE,recursive=FALSE)
goodFiles = c(1,2,4,8,11,13:15,17:19)

# find covariate data
varFiles = list.files(covarDir,pattern=".csv",full.names=TRUE,recursive=FALSE)
staticFiles = list.files(covarDir,pattern="Geo",full.names=TRUE)
staticFiles = c(staticFiles,list.files(covarDir,pattern="Slope",full.names=TRUE))

for (i in goodFiles) {
  
  # Initialize temp data frame
  thisSpecies = numeric()
  
  # Load in a file
  thisFile = data.frame(read.csv(presFiles[i]))
  colnames(thisFile) = c("Date",sites)
  
  # Get species from file name
  species = str_remove(presFiles[i],"_Daily.csv")
  species = str_remove(species,paste(presDir,'/',sep=""))
  if (str_detect(species,"Atl")){
    species = "Gervais"
  }
  
  # stack presence data from all sites
  thisSpecies = data.frame(Date=as.numeric(rep(as.numeric(as.Date(thisFile[,1],"%d-%b-%Y",origin="1970-01-01")),times=10)),
                           Pres=as.numeric(stack(thisFile[,2:11])[,1]),
                           Site=as.character(rep(sites[1:10],each=dim(thisFile)[1])))

  # add GS latitudinal position time series
  load(staticFiles[1])
  # Interpolate missing data
  # GSLat = approx(masterData.Time,masterData.Frontal,as.numeric(as.Date(thisFile$Date,"%d-%b-%Y",origin="1970-01-01")),method="linear")
  GSLat = approx(masterData.Time,masterData.Frontal,seq(min(masterData.Time),max(masterData.Time),by=1),method="linear")
  
  
  # load Dist to GS
  load(staticFiles[2])
  # load Slope & Aspect
  load(staticFiles[3])
  
  thisSpecies$GSLat = NA
  # for (j in 1:length(lags)){
  #   eval(parse(text=paste('thisSpecies$GSLatLag',lags[j],' = NA',sep="")))
  # }
  thisSpecies$GSDist = NA
  thisSpecies$Slope = NA
  thisSpecies$Aspect = NA
  
  for (n in 1:length(sites)){
    # interpolate GS Dist to fill in missing time stamps
    # GSdist = approx(masterData.Time,masterData.GeoDist[n,],as.numeric(as.Date(thisFile$Date,"%d-%b-%Y",origin="1970-01-01")),method="linear")
    GSdist = approx(masterData.Time,masterData.GeoDist[n,],seq(min(masterData.Time),max(masterData.Time),by=1),method="linear")
    
    # Plug data points in at correct time stamps
    siteInd = which(!is.na(str_match(thisSpecies$Site,sites[n])))
    putWhere = match(GSLat$x,thisSpecies$Date[siteInd])
    thisSpecies$GSLat[siteInd[putWhere[!is.na(putWhere)]]] = GSLat$y[-which(is.na(putWhere))]
    thisSpecies$GSDist[siteInd[putWhere[!is.na(putWhere)]]] = GSdist$y[-which(is.na(putWhere))]
    
    # # Create time lagged GSLat vectors
    # startInd = which(GSLat$x==as.Date('2016-05-01',origin="1970-01-01"))
    # fullLength = length(siteInd)
    #   for (k in lags){
    #     lagInd = startInd-k
    #     eval(parse(text=paste('thisSpecies$GSLatLag',as.character(k),'[siteInd] = GSLat$y[(lagInd):(lagInd+fullLength-1)]',sep="")))
    #   }
    
    if (n!=2 & n!=7){
    thisSpecies$Slope[siteInd] = slopeMat[n,1]
    thisSpecies$Aspect[siteInd] = aspectMat[n,1]
    } else if (n==2) {
      before = which(thisSpecies$Date[siteInd]<OC_change)
      thisSpecies$Slope[siteInd[before]] = slopeMat[n,1]
      thisSpecies$Slope[siteInd[-before]] = slopeMat[n,2]
      thisSpecies$Aspect[siteInd[before]] = aspectMat[n,1]
      thisSpecies$Aspect[siteInd[-before]] = aspectMat[n,2]
    } else if (n==7) {
      before = which(thisSpecies$Date[siteInd]<HAT_change)
      thisSpecies$Slope[siteInd[before]] = slopeMat[n,1]
      thisSpecies$Slope[siteInd[-before]] = slopeMat[n,2]
      thisSpecies$Aspect[siteInd[before]] = aspectMat[n,1]
      thisSpecies$Aspect[siteInd[-before]] = aspectMat[n,2]
    }
  }
  
  # Remove NA rows (no-effort days)
  noDat = which(thisSpecies$Pres=="NaN")
  thisSpecies = thisSpecies[-noDat,]
  rownames(thisSpecies) = seq(length=nrow(thisSpecies))
  
  if (i==11) {   # If Risso's, save data for later combination w UD36
    GgTemp = thisSpecies
    next
  } else if (i == 19) {  # If species is UD36, combine w Risso's
    thisSpecies = rbind(GgTemp,thisSpecies)
    # thisSpecies = thisSpecies %>% group_by(Date,Site) %>% summarise(Pres=sum(Pres),
    #                                                                 GSLat=mean(GSLat),
    #                                                                 GSLatLag7=mean(GSLatLag7),
    #                                                                 GSLatLag14=mean(GSLatLag14),
    #                                                                 GSLatLag21=mean(GSLatLag21),
    #                                                                 GSLatLag28=mean(GSLatLag28),
    #                                                                 GSLatLag42=mean(GSLatLag42),
    #                                                                 GSLatLag56=mean(GSLatLag56),
    #                                                                 GSDist=mean(GSDist),
    #                                                                 Slope=mean(Slope),
    #                                                                 Aspect=mean(Aspect))
    thisSpecies = thisSpecies %>% group_by(Date,Site) %>% summarise(Pres=sum(Pres),
                                                                    GSLat=mean(GSLat),
                                                                    GSDist=mean(GSDist),
                                                                    Slope=mean(Slope),
                                                                    Aspect=mean(Aspect))
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
        
        # for (m in 1:length(lags)){
          # # find covar cols with data for this site & depth & lag
          # thisSiteDepth = which(!is.na(str_match(colnames(thisVar),
          #                                        paste(sites[k],as.character(depths[l]),sep=""))) 
          #                       & str_detect(colnames(thisVar),paste('Lag',lags[m],sep="")))
          
        # find covar cols with data for this site & depth
        # thisSiteDepth = which(!is.na(str_match(colnames(thisVar),
        #                                        paste(sites[k],as.character(depths[l]),sep=""))))
        # 
        
          # add data to master data frame
          # if(k==1){
          # eval(parse(text=paste('thisSpecies$',abbrev,as.character(depths[l]),'Lag',lags[m],' = NA',sep="")))}
          # eval(parse(text=paste('thisSpecies$',abbrev,as.character(depths[l]),'Lag',
          #                       lags[m],'[siteInd[putWhere[!is.na(putWhere)]]] = thisVar[-which(is.na(putWhere)),thisSiteDepth]',sep="")))
        
        # add data to master data frame
        # if(k==1){
        #   eval(parse(text=paste('thisSpecies$',abbrev,as.character(depths[l]),' = NA',sep="")))}
        # eval(parse(text=paste('thisSpecies$',abbrev,as.character(depths[l]),
        #                       '[siteInd[putWhere[!is.na(putWhere)]]] = thisVar[-which(is.na(putWhere)),thisSiteDepth]',sep="")))
        
        # }
      }
    }
  }
  
  thisSpecies = apply(thisSpecies,2,as.character)
  saveName = paste(outDir,'/',species,'_masterDF.csv',sep="")
  write.csv(thisSpecies,saveName,row.names=FALSE)
  
}
