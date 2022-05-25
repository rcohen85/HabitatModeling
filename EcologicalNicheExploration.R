library(stringr)
library(dplyr)
library(lubridate)
library(sm)

inDir = 'J:/Chpt_3/ModelData'
outDir = 'J:/Chpt_3/EcologicalNichePlots'
fileList = list.files(path=inDir,pattern=paste('_masterDF.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)
covarList = c("EKE0","Temp0","Sal0","Temp200","Temp300","Temp400","Sal200","Sal300",
              "Sal400","Temp700","Sal700","SSH0","Chl0","FSLE0","VelMag0","VelMag200",
              "VelMag300","VelMag400","VelMag700","VelAsp0","VelAsp200","VelAsp300",
              "VelAsp400","VelAsp700")
stDt = as.Date("2016-05-01")
edDt = as.Date("2019-04-30")
allDates = stDt:edDt
weekDates = seq.Date(stDt,edDt,by=7)

sites = sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS')
OptImportanceStats = list()
# statNames = sort(paste(covarList,rep(c("_Opt","_Coef","_Pval","_AdjR","_Rho"),each=length(covarList)),sep=""))
statNames = sort(paste(covarList,rep(c("_Opt","_DistRho","_PresRho"),each=length(covarList)),sep=""))


for (i in 1:length(fileList)){
  
  spec = str_remove(str_remove(fileList[i],'_masterDF.csv'),paste(inDir,'/',sep=""))
  data = data.frame(read.csv(fileList[i]))
  data$Pres = round(data$Pres)
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(outDir,'/',spec,sep=""))){
    dir.create(paste(outDir,'/',spec,sep=""))
  }
  weeklyDF = as.numeric()

  # OptImportanceStats[[spec]]=data.frame(matrix(nrow=length(sites),ncol=length(covarList)*5))
  OptImportanceStats[[spec]]=data.frame(matrix(nrow=length(sites),ncol=length(covarList)*3))
  
  colnames(OptImportanceStats[[spec]]) = statNames
  
  for (l in 1:length(sites)) {
    
    # Create dataframe to hold data (or NAs) for all dates
    fullDatesDF = data.frame(matrix(nrow=length(allDates), ncol=43))
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
    
    summaryData$Site = sites[l]
    
    # find bins with presence (don't plot all the 0's)
    presWeeks = which(summaryData$Pres>0)
    
    
    for (j in 1:length(covarList)){  # for each covar of interest
      
      # calculate weekly average for this covar
      eval(parse(text=paste('thisCovar=fullDatesDF%>%group_by(WeekID)%>%summarize(',covarList[j],'=mean(',covarList[j],',na.rm=TRUE))',sep="")))
      eval(parse(text='summaryData[[covarList[j]]]=unlist(thisCovar[,2])'))
      
      if (i==1){
        
        png(paste(outDir,"/Full_",covarList[j],"Dist_at_",fullDatesDF$Site[1],".png",sep=""),width=350,height=350)
        plot(density(summaryData[[covarList[j]]],na.rm=TRUE),
             xlab=paste(covarList[j]),
             ylab='Density',
             main=paste(covarList[j]," at ",fullDatesDF$Site[1],sep=""))
        while (dev.cur()>1) {dev.off()}
        
      }
      
      if (sum(summaryData$Pres>0)>=5){ # only plot if sufficient presence at this site
        
        # fit line to pres vs. covar
        DatfitLine = lm(summaryData$Pres[presWeeks]~summaryData[[covarList[j]]][presWeeks])
        presCor = cor(summaryData$Pres[presWeeks],summaryData[[covarList[j]]][presWeeks])
        
        # calculate optimum value of covar as mean weighted by presence values
        optCovar = sum(summaryData$Pres[presWeeks]*summaryData[[covarList[j]]][presWeeks])/sum(summaryData$Pres[presWeeks])
        optDist = summaryData[[covarList[j]]][presWeeks]-optCovar
        
        # fit line to pres vs. dist from opt
        OptfitLine = lm(summaryData$Pres[presWeeks]~abs(optDist))
        
        # calculate correlation between presence and dist from opt
        distCor = cor(summaryData$Pres[presWeeks],abs(optDist))
        
        # add optimum value and importance stats to matrix
        OptImportanceStats[[spec]][[paste(covarList[j],"_Opt",sep="")]][l] = optCovar
        # OptImportanceStats[[spec]][[paste(covarList[j],"_Coef",sep="")]][l] = coef(OptfitLine)[2]
        # OptImportanceStats[[spec]][[paste(covarList[j],"_AdjR",sep="")]][l] = summary(OptfitLine)$adj.r.squared
        # OptImportanceStats[[spec]][[paste(covarList[j],"_Pval",sep="")]][l] = summary(OptfitLine)$coefficients[2,4]
        OptImportanceStats[[spec]][[paste(covarList[j],"_DistRho",sep="")]][l] = distCor
        OptImportanceStats[[spec]][[paste(covarList[j],"_PresRho",sep="")]][l] = presCor
        
        # Plot and save
        png(paste(outDir,"/",spec,"/Weekly",covarList[j],"_at_",fullDatesDF$Site[1],".png",sep=""),width=800,height=350)
        
        par(mfrow=c(1,3))
        # Plot presence vs. covar values w. line indicating value of "optimum"
        textX = (0.5*(max(summaryData[[covarList[j]]][presWeeks])-min(summaryData[[covarList[j]]][presWeeks])))+min(summaryData[[covarList[j]]][presWeeks])
        textY = 0.9*max(summaryData$Pres[presWeeks])
        plot(summaryData[[covarList[j]]][presWeeks],summaryData$Pres[presWeeks],type="p",
             ylab='Weekly Presence',
             xlab=paste(covarList[j]),
             cex.lab=1.5,
             cex.axis=1.5)
        # add fit line
        abline(DatfitLine,col="red")
        # Add line showing optimum value for this covar
        abline(v=optCovar,col="blue",lty=2)
        text(labels=paste("AdjR2 = ",round(summary(DatfitLine)$adj.r.squared,digits=4),"\n",
                          "P-value = ",round(summary(DatfitLine)$coefficients[2,4],digits=4),"\n",
                          "Rho = ",round(presCor,digits=4),sep=""),
             x=textX,
             y=textY,
             cex=1.5)
        
        # Plot presence vs. distance from optimum covar value
        textX = (0.5*(max(abs(optDist))-min(abs(optDist))))+min(abs(optDist))
        textY = 0.9*max(summaryData$Pres[presWeeks])
        plot(abs(optDist),summaryData$Pres[presWeeks],type="p",xlab='Distance from Optimum',
             ylab='Weekly Presence',cex.lab=1.5,cex.axis=1.5)
        #Add a fit line to this plot
        abline(OptfitLine,col="red")
        mtext(paste(spec,'and',covarList[j],'at',fullDatesDF$Site[1]),
              side = 3,
              line = - 2,
              outer = TRUE,
              cex=1.5)
        text(labels=paste("AdjR2 = ",round(summary(OptfitLine)$adj.r.squared,digits=4),"\n",
                          "P-value = ",round(summary(OptfitLine)$coefficients[2,4],digits=4),"\n",
                          "Rho = ",round(distCor,digits=4),sep=""),
             x=textX,
             y=textY,
             cex=1.5)
        
        # Plot distribution of presence covar values
        plot(density(summaryData[[covarList[j]]][presWeeks],na.rm=TRUE),
             xlab=paste(covarList[j]),
             ylab='Density',
             main=character(0),
             cex.lab=1.5,
             cex.axis=1.5)
        # Add line showing optimum value for this covar
        abline(v=optCovar,col="red",lty=2)
        
        while (dev.cur()>1) {dev.off()}
        
      }
    }
    
    if (sum(summaryData$Pres>0)>=5){
    weeklyDF = rbind(weeklyDF,summaryData)
    }
    
  }
  
  # Plot optimums for each covar across all sites
  for (j in 1:length(covarList)){
    
    presWeeks = which(weeklyDF$Pres>0)
    
    # calculate optimum value of covar as mean weighted by presence values
    optCovar = sum(weeklyDF$Pres[presWeeks]*weeklyDF[[covarList[j]]][presWeeks])/sum(weeklyDF$Pres[presWeeks])
    optDist = weeklyDF[[covarList[j]]][presWeeks]-optCovar
    
    # fit line to pres vs. dist from opt
    OptfitLine = lm(weeklyDF$Pres[presWeeks]~abs(optDist))
    
    # calculate correlation between presence and dist from opt
    distCor = cor(weeklyDF$Pres[presWeeks],abs(optDist))
    
    # Plot and save
    png(paste(outDir,"/",spec,"/Weekly",covarList[j],"_allSites.png",sep=""),width=800,height=350)
    
    par(mfrow=c(1,3))
    # Plot presence vs. covar values w. lines indicating value of "optimums"
    plot(weeklyDF[[covarList[j]]][presWeeks],weeklyDF$Pres[presWeeks],type="p",
         xlab=paste(covarList[j]),
         ylab='Weekly Presence',
         cex.lab=1.5,
         cex.axis=1.5)
    # Add lines showing site-specific optimum values
    abline(v=OptImportanceStats[[spec]][[paste(covarList[j],"_Opt",sep="")]],col="blue",lty=2)
    # Add line showing optimum value for this covar across sites
    abline(v=optCovar,col="red",lty=2)
    
    # Plot presence vs. distance from optimum covar value
    textX = (0.5*(max(abs(optDist))-min(abs(optDist))))+min(abs(optDist))
    textY = 0.9*max(weeklyDF$Pres[presWeeks])
    plot(abs(optDist),weeklyDF$Pres[presWeeks],type="p",xlab='Distance from Optimum',ylab='Weekly Presence',
         cex.lab=1.5,
         cex.axis=1.5)
    #Add a fit line to this plot
    abline(OptfitLine,col="red")
    mtext(paste(spec,'and',covarList[j],'at all sites'),
          side = 3,
          line = - 2,
          outer = TRUE,
          cex=1.5)
    text(labels=paste("AdjR2 = ",round(summary(OptfitLine)$adj.r.squared,digits=4),"\n",
                      "P-value = ",round(summary(OptfitLine)$coefficients[2,4],digits=4),"\n",
                      "Rho = ",round(distCor,digits=4),sep=""),
         x=textX,
         y=textY,
         cex=1.5)
    
    # Plot distribution of presence covar values
    # sm.density.compare(weeklyDF[[covarList[j]]][presWeeks],weeklyDF$Site[presWeeks],
    #                         xlab=paste(covarList[j]),
    #                         ylab='Density',
    #                         main=character(0),
    #                         cex.lab=1.5,
    #                         cex.axis=1.5)
    plot(density(weeklyDF[[covarList[j]]][presWeeks],na.rm=TRUE),
         xlab=paste(covarList[j]),
         ylab='Density',
         main=character(0),
         cex.lab=1.5,
         cex.axis=1.5)
    # Add lines showing site-specific optimum values
    abline(v=OptImportanceStats[[spec]][[paste(covarList[j],"_Opt",sep="")]],col="blue",lty=2)
    # Add line showing optimum value for this covar across sites
    abline(v=optCovar,col="red",lty=2)
    
    while (dev.cur()>1) {dev.off()}
    
    
  }
  
}

save(OptImportanceStats,file=paste(outDir,'OptimumImportance'))
