library(stringr)
library(dplyr)
library(lubridate)
inDir = 'E:/ModelingCovarData/Master_DFs'
outDir = 'E:/ModelingCovarData/EcologicalNichePlots'
fileList = list.files(path=inDir,pattern=paste('_masterDF.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)
covarList = c("Temp0","Sal0","Temp300","Sal300","Temp700","Sal700","SSH0","Chl0","FSLE0")
stDt = as.Date("2016-05-01")
edDt = as.Date("2019-04-30")
allDates = stDt:edDt
weekDates = seq.Date(stDt,edDt,by=7)

sites = sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS')

for (i in 1:length(fileList)){
  
  spec = str_remove(str_remove(fileList[i],'_masterDF.csv'),paste(inDir,'/',sep=""))
  data = data.frame(read.csv(fileList[i]))
  data$Pres = round(data$Pres)
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(outDir,'/',spec,sep=""))){
    dir.create(paste(outDir,'/',spec,sep=""))
  }
  weeklyDF = as.numeric()
  
  for (l in 1:length(sites)) {
    
    # Create dataframe to hold data (or NAs) for all dates
    fullDatesDF = data.frame(matrix(nrow=length(allDates), ncol=258))
    fullDatesDF[,1] = allDates
    
    # sort the observations we have for this site into the full date range
    thisSite = which(!is.na(str_match(data$Site,sites[l])))
    matchRow = match(data$Date[thisSite],allDates)
    fullDatesDF[matchRow,2:dim(data)[2]] = data[thisSite,2:dim(data)[2]]
    
    colnames(fullDatesDF) = colnames(data)
    
    # create grouping variable
    weekID = rep(1:length(weekDates),each=7)
    weekID = weekID[1:dim(fullDatesDF)[1]]
    fullDatesDF$WeekID = weekID
    
    summaryData = fullDatesDF %>%
      group_by(WeekID) %>%
      summarize(Pres=sum(Pres,na.rm=TRUE))
    
    if (sum(summaryData$Pres>0)>=5){ # only plot if sufficient presence at this site
      for (j in 1:length(covarList)){  # for each covar of interest
        
        # calculate weekly average for this covar
        eval(parse(text=paste('thisCovar=fullDatesDF%>%group_by(WeekID)%>%summarize(',covarList[j],'=mean(',covarList[j],',na.rm=TRUE))',sep="")))
        eval(parse(text='summaryData[[covarList[j]]]=unlist(thisCovar[,2])'))
        
        # TO DO: find bins with presence (don't plot all the 0's)
        presWeeks = which(summaryData$Pres>0)
        
        # calculate optimum value of covar as mean weighted by presence values
        optCovar = sum(summaryData$Pres[presWeeks]*summaryData[[covarList[j]]][presWeeks])/sum(summaryData$Pres[presWeeks])
        optDist = summaryData[[covarList[j]]][presWeeks]-optCovar
        
        # Plot and save
        png(paste(outDir,"/",spec,"/Weekly",covarList[j],"_at_",fullDatesDF$Site[1],".png",sep=""),width=800,height=400)
        
        par(mfrow=c(1,2))
        plot(summaryData[[covarList[j]]][presWeeks],summaryData$Pres[presWeeks],type="p",
             ylab='Weekly Presence',
             xlab=paste(covarList[j]))
        # Add line showing optimum value for this covar
        abline(v=optCovar,col="red",lty=2)
        
        plot(abs(optDist),summaryData$Pres[presWeeks],type="p",xlab='Distance from Optimum',ylab='Weekly Presence')
        #Add a fit line to this plot
        abline(lm(summaryData$Pres[presWeeks]~abs(optDist)),col="red")
        mtext(paste(spec,'and',covarList[j],'at',fullDatesDF$Site[1]),
              side = 3,
              line = - 2,
              outer = TRUE,
              cex=1.5)
        
        while (dev.cur()>1) {dev.off()}
        
      }
      
      weeklyDF = rbind(weeklyDF,summaryData)
    }
  }
  
  # Plot optimums for each covar across all sites
  for (j in 1:length(covarList)){
    
    presWeeks = which(weeklyDF$Pres>0)
    
    # calculate optimum value of covar as mean weighted by presence values
    optCovar = sum(weeklyDF$Pres[presWeeks]*weeklyDF[[covarList[j]]][presWeeks])/sum(weeklyDF$Pres[presWeeks])
    optDist = weeklyDF[[covarList[j]]][presWeeks]-optCovar
    
    # Plot and save
    png(paste(outDir,"/",spec,"/Weekly",covarList[j],"_allSites.png",sep=""),width=800,height=400)
    
    par(mfrow=c(1,2))
    plot(weeklyDF[[covarList[j]]][presWeeks],weeklyDF$Pres[presWeeks],type="p",
         xlab=paste(covarList[j]),
         ylab='Weekly Presence')
    # Add line showing optimum value for this covar
    abline(v=optCovar,col="red",lty=2)
    
    plot(abs(optDist),weeklyDF$Pres[presWeeks],type="p",xlab='Distance from Optimum',ylab='Weekly Presence')
    #Add a fit line to this plot
    abline(lm(weeklyDF$Pres[presWeeks]~abs(optDist)),col="red")
    mtext(paste(spec,'and',covarList[j],'at all sites'),
          side = 3,
           line = - 2,
           outer = TRUE,
          cex=1.5)
    
    while (dev.cur()>1) {dev.off()}
    
    
  }
  
}
