library(stringr)
library(dplyr)
library(lubridate)
library(sm)
library(pracma)
library(ggplot2)
library(gridExtra)
library(gridGraphics)
library(rstatix)
library(expss)
library(corrplot)
library(cluster)
library(ggbiplot)

inDir = 'J:/Chpt_3/ModelData'
outDir = 'J:/Chpt_3/EcologicalNichePlots'
# fileList = list.files(path=inDir,pattern=paste('_masterDF.csv',sep=""),
#                       full.names=TRUE,recursive=FALSE,
#                       include.dirs=FALSE,no..=TRUE)
fileList = list.files(path=getwd(),pattern=paste('_masterDF.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)


## Plot presence/covar correlations, "optimum" values, etc for each species --------------
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS')
covarList = c("EKE0","Temp0","Sal0","Temp200","Temp300","Temp400","Sal200","Sal300",
              "Sal400","Temp700","Sal700","SSH0","Chl0","FSLE0","VelMag0","VelMag200",
              "VelMag300","VelMag400","VelMag700","VelAsp0","VelAsp200","VelAsp300",
              "VelAsp400","VelAsp700","AEddyDist0","CEddyDist0","GSLat","GSDist")
stDt = as.Date("2016-05-01")
edDt = as.Date("2019-04-30")
allDates = stDt:edDt
weekDates = seq.Date(stDt,edDt,by=7)
WeekOptStats = list()
DayOptStats = list()
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
  
  WeekOptStats[[spec]]=data.frame(matrix(nrow=length(sites),ncol=length(covarList)*3))
  DayOptStats[[spec]]=data.frame(matrix(nrow=length(sites),ncol=length(covarList)*3))
  
  colnames(WeekOptStats[[spec]]) = statNames
  colnames(DayOptStats[[spec]]) = statNames
  
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
    
    # sum presence in each week
    summaryData = fullDatesDF %>%
      group_by(WeekID) %>%
      summarize(Pres=sum(Pres,na.rm=TRUE))
    
    # normalize by effort
    effDF = data.frame(count=rep(1,length(allDates)),weekID=weekID)
    effDF$count[which(is.na(fullDatesDF$Pres))] = 0
    propEff = effDF %>%
      group_by(weekID) %>%
      summarize(sum(count))
    summaryData$propEff = unlist(propEff[,2])/7
    summaryData$Pres = summaryData$Pres*(1/summaryData$propEff)
    
    summaryData$Site = sites[l]
    
    # find bins with presence (don't plot all the 0's)
    presWeeks = which(summaryData$Pres>0)
    
    
    for (j in 1:length(covarList)){  # for each covar of interest
      
      # calculate weekly average for this covar
      eval(parse(text=paste('thisCovar=fullDatesDF%>%group_by(WeekID)%>%summarize(',covarList[j],'=mean(',covarList[j],',na.rm=TRUE))',sep="")))
      eval(parse(text='summaryData[[covarList[j]]]=unlist(thisCovar[,2])'))
      
      if (i==1){
        
        png(paste(outDir,"/FullWeekly_",covarList[j],"Dist_at_",sites[l],".png",sep=""),width=350,height=350)
        plot(density(summaryData[[covarList[j]]],na.rm=TRUE),
             xlab=paste(covarList[j]),
             ylab='Density',
             main=paste(covarList[j]," at ",sites[l],sep=""))
        while (dev.cur()>1) {dev.off()}
        
        png(paste(outDir,"/FullDaily_",covarList[j],"Dist_at_",sites[l],".png",sep=""),width=350,height=350)
        plot(density(data[[covarList[j]]],na.rm=TRUE),
             xlab=paste(covarList[j]),
             ylab='Density',
             main=paste(covarList[j]," at ",sites[l],sep=""))
        while (dev.cur()>1) {dev.off()}
        
      }
      
      if (sum(summaryData$Pres>0,na.rm=TRUE)>=5){ # only plot if sufficient presence at this site
        
        # fit line to pres vs. covar
        weekDatfitLine = lm(summaryData$Pres[presWeeks]~summaryData[[covarList[j]]][presWeeks])
        weekpresCor = cor(summaryData$Pres[presWeeks],summaryData[[covarList[j]]][presWeeks])
        
        # calculate optimum value of covar as mean weighted by presence values
        weekoptCovar = sum(summaryData$Pres[presWeeks]*summaryData[[covarList[j]]][presWeeks])/sum(summaryData$Pres[presWeeks])
        weekoptDist = summaryData[[covarList[j]]][presWeeks]-weekoptCovar
        
        # fit line to pres vs. dist from opt
        weekOptfitLine = lm(summaryData$Pres[presWeeks]~abs(weekoptDist))
        
        # calculate correlation between presence and dist from opt
        weekdistCor = cor(summaryData$Pres[presWeeks],abs(weekoptDist))
        
        # add optimum value and importance stats to matrix
        WeekOptStats[[spec]][[paste(covarList[j],"_Opt",sep="")]][l] = weekoptCovar
        WeekOptStats[[spec]][[paste(covarList[j],"_DistRho",sep="")]][l] = weekdistCor
        WeekOptStats[[spec]][[paste(covarList[j],"_PresRho",sep="")]][l] = weekpresCor
        
        # Plot and save weekly optimums
        png(paste(outDir,"/",spec,"/Weekly",covarList[j],"_at_",sites[l],".png",sep=""),width=800,height=350)
        
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
        abline(weekDatfitLine,col="red")
        # Add line showing optimum value for this covar
        abline(v=weekoptCovar,col="blue",lty=2)
        text(labels=paste("AdjR2 = ",round(summary(weekDatfitLine)$adj.r.squared,digits=4),"\n",
                          "P-value = ",round(summary(weekDatfitLine)$coefficients[2,4],digits=4),"\n",
                          "Rho = ",round(weekpresCor,digits=4),sep=""),
             x=textX,
             y=textY,
             cex=1.5)
        
        # Plot presence vs. distance from optimum covar value
        textX = (0.5*(max(abs(weekoptDist))-min(abs(weekoptDist))))+min(abs(weekoptDist))
        textY = 0.9*max(summaryData$Pres[presWeeks])
        plot(abs(weekoptDist),summaryData$Pres[presWeeks],type="p",xlab='Distance from Optimum',
             ylab='Weekly Presence',cex.lab=1.5,cex.axis=1.5)
        #Add a fit line to this plot
        abline(weekOptfitLine,col="red")
        mtext(paste(spec,'and',covarList[j],'at',sites[l]),
              side = 3,
              line = - 2,
              outer = TRUE,
              cex=1.5)
        text(labels=paste("AdjR2 = ",round(summary(weekOptfitLine)$adj.r.squared,digits=4),"\n",
                          "P-value = ",round(summary(weekOptfitLine)$coefficients[2,4],digits=4),"\n",
                          "Rho = ",round(weekdistCor,digits=4),sep=""),
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
        abline(v=weekoptCovar,col="red",lty=2)
        
        while (dev.cur()>1) {dev.off()}
        
      }
      
      if (sum(data$Pres[thisSite]>0,na.rm=TRUE)>=25){
        
        presDays = which(data$Pres[thisSite]>0)
        
        # fit line to pres vs. covar
        dayDatfitLine = lm(data$Pres[thisSite[presDays]]~data[[covarList[j]]][thisSite[presDays]])
        daypresCor = cor(data$Pres[thisSite[presDays]],data[[covarList[j]]][thisSite[presDays]])
        
        # calculate optimum value of covar as mean weighted by presence values
        dayoptCovar = sum(data$Pres[thisSite[presDays]]*data[[covarList[j]]][thisSite[presDays]])/sum(data$Pres[thisSite[presDays]])
        dayoptDist = data[[covarList[j]]][thisSite[presDays]]-dayoptCovar
        
        # fit line to pres vs. dist from opt
        dayOptfitLine = lm(data$Pres[thisSite[presDays]]~abs(dayoptDist))
        
        # calculate correlation between presence and dist from opt
        daydistCor = cor(data$Pres[thisSite[presDays]],abs(dayoptDist))
        
        # add optimum value and importance stats to matrix
        DayOptStats[[spec]][[paste(covarList[j],"_Opt",sep="")]][l] = dayoptCovar
        DayOptStats[[spec]][[paste(covarList[j],"_DistRho",sep="")]][l] = daydistCor
        DayOptStats[[spec]][[paste(covarList[j],"_PresRho",sep="")]][l] = daypresCor
        
        # Plot and save daily optimums
        png(paste(outDir,"/",spec,"/Daily",covarList[j],"_at_",sites[l],".png",sep=""),width=800,height=350)
        
        par(mfrow=c(1,3))
        # Plot presence vs. covar values w. line indicating value of "optimum"
        textX = (0.5*(max(data[[covarList[j]]][thisSite[presDays]])-min(data[[covarList[j]]][thisSite[presDays]])))+min(data[[covarList[j]]][thisSite[presDays]])
        textY = 0.9*max(data$Pres[thisSite[presDays]])
        plot(data[[covarList[j]]][thisSite[presDays]],data$Pres[thisSite[presDays]],type="p",
             ylab='Daily Presence',
             xlab=paste(covarList[j]),
             cex.lab=1.5,
             cex.axis=1.5)
        # add fit line
        abline(dayDatfitLine,col="red")
        # Add line showing optimum value for this covar
        abline(v=dayoptCovar,col="blue",lty=2)
        text(labels=paste("AdjR2 = ",round(summary(dayDatfitLine)$adj.r.squared,digits=4),"\n",
                          "P-value = ",round(summary(dayDatfitLine)$coefficients[2,4],digits=4),"\n",
                          "Rho = ",round(daypresCor,digits=4),sep=""),
             x=textX,
             y=textY,
             cex=1.5)
        
        # Plot presence vs. distance from optimum covar value
        textX = (0.5*(max(abs(dayoptDist))-min(abs(dayoptDist))))+min(abs(dayoptDist))
        textY = 0.9*max(data$Pres[thisSite[presDays]])
        plot(abs(dayoptDist),data$Pres[thisSite[presDays]],type="p",xlab='Distance from Optimum',
             ylab='Daily Presence',cex.lab=1.5,cex.axis=1.5)
        #Add a fit line to this plot
        abline(dayOptfitLine,col="red")
        mtext(paste(spec,'and',covarList[j],'at',sites[l]),
              side = 3,
              line = - 2,
              outer = TRUE,
              cex=1.5)
        text(labels=paste("AdjR2 = ",round(summary(dayOptfitLine)$adj.r.squared,digits=4),"\n",
                          "P-value = ",round(summary(dayOptfitLine)$coefficients[2,4],digits=4),"\n",
                          "Rho = ",round(daydistCor,digits=4),sep=""),
             x=textX,
             y=textY,
             cex=1.5)
        
        # Plot distribution of presence covar values
        plot(density(data[[covarList[j]]][thisSite[presDays]],na.rm=TRUE),
             xlab=paste(covarList[j]),
             ylab='Density',
             main=character(0),
             cex.lab=1.5,
             cex.axis=1.5)
        # Add line showing optimum value for this covar
        abline(v=dayoptCovar,col="red",lty=2)
        
        while (dev.cur()>1) {dev.off()}
        
      }
    }
    
    if (sum(summaryData$Pres>0,na.rm=TRUE)>=5){
      weeklyDF = rbind(weeklyDF,summaryData)
    }
    
  }
  
  # Plot daily and weekly optimums for each covar across all sites
  for (j in 1:length(covarList)){
    
    presWeeks = which(weeklyDF$Pres>0)
    
    # calculate correlation between presence and covar
    weekCor = cor(weeklyDF$Pres[presWeeks],weeklyDF[[covarList[j]]][presWeeks])
    weekFitLine = lm(weeklyDF$Pres[presWeeks]~weeklyDF[[covarList[j]]][presWeeks])
    
    # calculate optimum value of covar as mean weighted by presence values
    weekoptCovar = sum(weeklyDF$Pres[presWeeks]*weeklyDF[[covarList[j]]][presWeeks])/sum(weeklyDF$Pres[presWeeks])
    weekoptDist = weeklyDF[[covarList[j]]][presWeeks]-weekoptCovar
    
    # fit line to pres vs. dist from opt
    weekOptfitLine = lm(weeklyDF$Pres[presWeeks]~abs(weekoptDist))
    
    # calculate correlation between presence and dist from opt
    weekdistCor = cor(weeklyDF$Pres[presWeeks],abs(weekoptDist))
    
    # Plot and save
    png(paste(outDir,"/",spec,"/Weekly",covarList[j],"_allSites.png",sep=""),width=800,height=350)
    
    par(mfrow=c(1,3))
    # Plot presence vs. covar values w. lines indicating value of "optimums" and overall correlation
    textX = (0.5*(max(weeklyDF[[covarList[j]]][presWeeks])-min(weeklyDF[[covarList[j]]][presWeeks])))+min(weeklyDF[[covarList[j]]][presWeeks])
    textY = 0.9*max(weeklyDF$Pres[presWeeks])
    plot(weeklyDF[[covarList[j]]][presWeeks],weeklyDF$Pres[presWeeks],type="p",
         xlab=paste(covarList[j]),
         ylab='Weekly Presence',
         cex.lab=1.5,
         cex.axis=1.5)
    # Add lines showing site-specific optimum values
    abline(v=WeekOptStats[[spec]][[paste(covarList[j],"_Opt",sep="")]],col="blue",lty=2)
    # Add line showing optimum value for this covar across sites
    abline(v=weekoptCovar,col="red",lty=2)
    # add line showing overall relationship between presence & covar
    abline(weekFitLine,col="green")
    text(labels=paste("Rho = ",round(weekCor,digits=4),sep=""),
         x=textX,
         y=textY,
         cex=1.5)
    
    # Plot presence vs. distance from optimum covar value
    textX = (0.5*(max(abs(weekoptDist))-min(abs(weekoptDist))))+min(abs(weekoptDist))
    textY = 0.9*max(weeklyDF$Pres[presWeeks])
    plot(abs(weekoptDist),weeklyDF$Pres[presWeeks],type="p",xlab='Distance from Optimum',ylab='Weekly Presence',
         cex.lab=1.5,
         cex.axis=1.5)
    #Add a fit line to this plot
    abline(weekOptfitLine,col="red")
    mtext(paste(spec,'and',covarList[j],'at all sites'),
          side = 3,
          line = - 2,
          outer = TRUE,
          cex=1.5)
    text(labels=paste("AdjR2 = ",round(summary(weekOptfitLine)$adj.r.squared,digits=4),"\n",
                      "P-value = ",round(summary(weekOptfitLine)$coefficients[2,4],digits=4),"\n",
                      "Rho = ",round(weekdistCor,digits=4),sep=""),
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
    abline(v=WeekOptStats[[spec]][[paste(covarList[j],"_Opt",sep="")]],col="blue",lty=2)
    # Add line showing optimum value for this covar across sites
    abline(v=weekoptCovar,col="red",lty=2)
    
    while (dev.cur()>1) {dev.off()}
    
    presDays = which(data$Pres>0)
    
    # calculate correlation between presence and covar
    dayCor = cor(data$Pres[presDays],data[[covarList[j]]][presDays])
    dayFitLine = lm(data$Pres[presDays]~data[[covarList[j]]][presDays])
    
    # calculate optimum value of covar as mean weighted by presence values
    dayoptCovar = sum(data$Pres[presDays]*data[[covarList[j]]][presDays])/sum(data$Pres[presDays])
    dayoptDist = data[[covarList[j]]][presDays]-dayoptCovar
    
    # fit line to pres vs. dist from opt
    dayOptfitLine = lm(data$Pres[presDays]~abs(dayoptDist))
    
    # calculate correlation between presence and dist from opt
    daydistCor = cor(data$Pres[presDays],abs(dayoptDist))
    
    # Plot and save
    png(paste(outDir,"/",spec,"/Daily",covarList[j],"_allSites.png",sep=""),width=800,height=350)
    
    par(mfrow=c(1,3))
    # Plot presence vs. covar values w. lines indicating value of "optimums"
    textX = (0.5*(max(data[[covarList[j]]][presDays])-min(data[[covarList[j]]][presDays])))+min(data[[covarList[j]]][presDays])
    textY = 0.9*max(data$Pres[presDays])
    plot(data[[covarList[j]]][presDays],data$Pres[presDays],type="p",
         xlab=paste(covarList[j]),
         ylab='Daily Presence',
         cex.lab=1.5,
         cex.axis=1.5)
    # Add lines showing site-specific optimum values
    abline(v=DayOptStats[[spec]][[paste(covarList[j],"_Opt",sep="")]],col="blue",lty=2)
    # Add line showing optimum value for this covar across sites
    abline(v=dayoptCovar,col="red",lty=2)
    # add line showing overall relationship between presence & covar
    abline(dayFitLine,col="green")
    text(labels=paste("Rho = ",round(dayCor,digits=4),sep=""),
         x=textX,
         y=textY,
         cex=1.5)
    
    # Plot presence vs. distance from optimum covar value
    textX = (0.5*(max(abs(dayoptDist))-min(abs(dayoptDist))))+min(abs(dayoptDist))
    textY = 0.9*max(data$Pres[presDays])
    plot(abs(dayoptDist),data$Pres[presDays],type="p",xlab='Distance from Optimum',ylab='Daily Presence',
         cex.lab=1.5,
         cex.axis=1.5)
    #Add a fit line to this plot
    abline(dayOptfitLine,col="red")
    mtext(paste(spec,'and',covarList[j],'at all sites'),
          side = 3,
          line = - 2,
          outer = TRUE,
          cex=1.5)
    text(labels=paste("AdjR2 = ",round(summary(dayOptfitLine)$adj.r.squared,digits=4),"\n",
                      "P-value = ",round(summary(dayOptfitLine)$coefficients[2,4],digits=4),"\n",
                      "Rho = ",round(daydistCor,digits=4),sep=""),
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
    plot(density(data[[covarList[j]]][presDays],na.rm=TRUE),
         xlab=paste(covarList[j]),
         ylab='Density',
         main=character(0),
         cex.lab=1.5,
         cex.axis=1.5)
    # Add lines showing site-specific optimum values
    abline(v=DayOptStats[[spec]][[paste(covarList[j],"_Opt",sep="")]],col="blue",lty=2)
    # Add line showing optimum value for this covar across sites
    abline(v=dayoptCovar,col="red",lty=2)
    
    while (dev.cur()>1) {dev.off()}
    
    
  }
  
}

save(WeekOptStats,DayOptStats,file=paste(outDir,'/','OptimumImportance.Rdata',sep=""))

## Compare species covar selections at specific sites ------------
# specs = c("UD26","Risso")
# sites = c("HZ","BC","NFC","HAT")
# covars = c("Chl0","FSLE0","SSH0","Sal0","Sal400","Temp0","Temp400","VelAsp0","VelAsp400","VelMag0","VelMag400","AEddyDist0","CEddyDist0")
# specs = c("Blainville","Gervais")
# sites = c("BS")
specs = c("Sowerby","Cuvier")
sites = c("HZ","BC","WC")
covars = c("Chl0","FSLE0","SSH0","Sal0","Sal700","Temp0","Temp700","VelAsp0","VelAsp700","VelMag0","VelMag700","AEddyDist0","CEddyDist0")
plotColors = c('#D62A1C','#1C3FD6')

files = c(str_which(fileList,specs[1]),str_which(fileList,specs[2]))

# if it doesn't already exist, create directory to save figures
if (!dir.exists(paste(outDir,'/',specs[1],"_",specs[2],sep=""))){
  dir.create(paste(outDir,'/',specs[1],"_",specs[2],sep=""))
}

data1 = data.frame(read.csv(fileList[files[1]]))
data1$Pres = round(data1$Pres)
data2 = data.frame(read.csv(fileList[files[2]]))
data2$Pres = round(data2$Pres)

for (j in 1:length(sites)){
  for (i in 1:length(covars)){

    # find data at this site
    siteInd1 = str_which(data1$Site,sites[j])
    siteInd2 = str_which(data2$Site,sites[j])

    # Divide covar range into bins
    varMin = min(min(data1[[covars[i]]][siteInd1]),min(data2[[covars[i]]][siteInd2]))
    varMin = varMin-abs(0.01*varMin)
    varMax = max(max(data1[[covars[i]]][siteInd1]),max(data2[[covars[i]]][siteInd2]))
    varMax = varMax+abs(0.01*varMax)
    varBins = seq(varMin,varMax,length.out=21)

    # Identify which bin each covar observation falls into
    binDat1 = histc(data1[[covars[i]]][siteInd1],varBins)
    binDat2 = histc(data2[[covars[i]]][siteInd2],varBins)

    # Calculate average and SD of presence values corresponding to covar obs in each bin
    presChar1 = data.frame(pres=data1$Pres[siteInd1],bin=binDat1$bin) %>%
      group_by(bin) %>%
      summarize(MeanPres=mean(pres,na.rm=TRUE),
                SD=sd(pres,na.rm=TRUE))
    presChar2 = data.frame(pres=data2$Pres[siteInd2],bin=binDat2$bin) %>%
      group_by(bin) %>%
      summarize(MeanPres=mean(pres,na.rm=TRUE),
                SD=sd(pres,na.rm=TRUE))

    # Normalize by max mean presence per species
    presChar1$NormMean = presChar1$MeanPres-min(presChar1$MeanPres)
    presChar1$NormMean = presChar1$NormMean/max(presChar1$NormMean)
    presChar1$NormSD = presChar1$SD-min(presChar1$MeanPres)
    presChar1$NormSD = presChar1$NormSD/max(presChar1$NormMean)
    presChar2$NormMean = presChar2$MeanPres-min(presChar2$MeanPres)
    presChar2$NormMean = presChar2$NormMean/max(presChar2$NormMean)
    presChar2$NormSD = presChar2$SD-min(presChar2$MeanPres)
    presChar2$NormSD = presChar2$NormSD/max(presChar2$NormMean)

    presChar1$binCenter = varBins[presChar1$bin]+0.5*diff(varBins[1:2])
    presChar2$binCenter = varBins[presChar2$bin]+0.5*diff(varBins[1:2])

    # Plot as overlaid bar charts
    barPlot = ggplot(
    )+geom_col(data=presChar1,
               aes(x=binCenter,y=NormMean),
               color=plotColors[1],
               fill=plotColors[1],
               alpha=0.5
    )+geom_col(data=presChar2,
               aes(x=binCenter,y=NormMean),
               color=plotColors[2],
               fill=plotColors[2],
               alpha=0.5
    )+coord_cartesian(xlim=c(varMin,varMax)
    )+labs(x=covars[i],y="Normalized Mean Presence",title=covars[i]
    )+theme_minimal()

    # make legend
    legDF = data.frame(x1=c(1,1),
                       x2=c(2,2),
                       y1=c(1,2),
                       y2=c(1.75,2.75),
                       cols=plotColors)
    legPlot = ggplot(legDF
    )+geom_rect(aes(ymin=y1,
                    ymax=y2,
                    xmin=x1,
                    xmax=x2,
                    color=cols,
                    fill=cols),
                alpha=0.5
    )+scale_fill_identity(
    )+scale_color_identity(
    )+geom_text(aes(x=x2+0.95,
                    y=y1+0.5,
                    label=specs),
                size=2.75
    )+guides(fill='none',color='none'
    )+coord_cartesian(xlim=c(1,4)
    )+theme_void()

    # png(file=paste(outDir,'/',specs[1],"_",specs[2],"/",covars[i],"_at_",sites[j],".png",sep=""),width = 500, height = 400, units = "px",res=125)
    png(file=paste(getwd(),'/EcologicalNichePlots/',specs[1],"_",specs[2],"_",covars[i],"_at_",sites[j],".png",sep=""),width = 500, height = 400, units = "px",res=125)
    grid.arrange(barPlot,legPlot,ncol=5,nrow=4,layout_matrix=rbind(c(rep(1,4),NA),
                                                               c(rep(1,4),2),
                                                               c(rep(1,4),NA),
                                                               c(rep(1,4),NA)))
    while (dev.cur()>1) {dev.off()}

  }
}

## Compare presence conditions across species using violinplots ------------------------

fileList = dir('J:/Chpt_3/ModelData',pattern=".csv")
covars = c("Chl0","EKE0","FSLE0","SSH0","Temp0","Temp200","Temp400","Temp700","Sal0","Sal200","Sal400","Sal700",
           "VelMag0","VelAsp0","CEddyDist0","AEddyDist0")
# specs = cbind(c("Blainville","Cuvier","Gervais","Kogia","Risso","Sowerby","SpermWhale","True","UD26","UD28"),
#               c("Md","Zc","Me","Kg","Gg","Mb","Pm","Mm","Gm","Dd"))
specs = cbind(rev(c("UD28","Risso","UD26","Blainville","Gervais","Cuvier","Sowerby","True","Kogia","SpermWhale")),
              rev(c("Dd","Gg","Gm","Md","Me","Zc","Mb","Mm","Kg","Pm")))

masterDF = list()
allNames = list()
for (i in 1:length(fileList)){
  
  data = data.frame(read.csv(fileList[i]))
  thisSpec = str_remove(fileList[i],'_masterDF.csv')
  specAbbrev = specs[which(str_detect(specs[,1],thisSpec)),2]
  
  # Round presence to get Poisson dist
  data$Pres = round(data$Pres)
  
  # create weekly time series to reduce autocorrelation
  stDt = as.Date("2016-05-01")
  edDt = as.Date("2019-04-30")
  sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS')
  allDates = stDt:edDt
  weekDates = seq.Date(stDt,edDt,by=7)
  weeklyDF = as.numeric()
  
  for (l in 1:length(sites)) {
    
    # Create dataframe to hold data (or NAs) for all dates
    fullDatesDF = data.frame(matrix(nrow=length(allDates), ncol=44))
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
    
    # sum presence in each week
    summaryData = fullDatesDF %>%
      group_by(WeekID) %>%
      summarize(Pres=sum(Pres,na.rm=TRUE))
    
    # normalize by effort
    effDF = data.frame(count=rep(1,length(allDates)),weekID=weekID)
    effDF$count[which(is.na(fullDatesDF$Pres))] = 0
    propEff = effDF %>%
      group_by(weekID) %>%
      summarize(sum(count))
    summaryData$propEff = unlist(propEff[,2])/7
    summaryData$Pres = summaryData$Pres*(1/summaryData$propEff)
    
    summaryData$Site = sites[l]
    summaryData$WeekDate = weekDates
    
    for (j in 4:(dim(data)[2])){
      var = names(data)[j]
      # calculate weekly average for this covar
      eval(parse(text=paste('thisCovar=fullDatesDF%>%group_by(WeekID)%>%summarize(',var,'=mean(',var,',na.rm=TRUE))',sep="")))
      eval(parse(text='summaryData[[var]]=unlist(thisCovar[,2])'))
      
    }
    weeklyDF = rbind(weeklyDF,summaryData)
  }
  
  
  if (i==1){
    masterDF$WeekDate = weeklyDF$WeekDate
    masterDF$Site = weeklyDF$Site
    allNames = c(allNames,"WeekDate","Site")
    for (j in 1:length(covars)){
      masterDF[[covars[j]]] = weeklyDF[[covars[j]]]
      allNames = c(allNames,covars[j])
    }
  }
  
  masterDF[[specAbbrev]] = weeklyDF$Pres
  allNames = c(allNames,specAbbrev)
}

masterDF = data.frame(masterDF)
write.csv(masterDF,'MasterWeeklyDF.csv',row.names=FALSE)

# make violin plots for each covar

masterDF = data.frame(read.csv('MasterWeeklyDF.csv'))
for (i in 1:length(covars)){
  plotDF = list()
  
  for (j in 1:dim(specs)[1]){
    presVals = as.numeric(masterDF[[covars[i]]][(masterDF[[specs[j,2]]]>0 & !is.na(masterDF[[specs[j,2]]]))])
    # plotDF = rbind(plotDF,cbind(presVals,rep(j,length(presVals))))
    plotDF[[specs[j,2]]] = presVals
  }
  
  plotDF = data.frame(stack(plotDF))
  
  xmin = min(masterDF[[covars[i]]],na.rm=TRUE)
  xmax = max(masterDF[[covars[i]]],na.rm=TRUE)
  
  ggplot(plotDF,
         # aes(x=ind,y=values,color=ind), # different color for each species
         aes(y=ind,x=values), # consistent color across species
  )+geom_violin(color="#43a5f0",
                fill="#43a5f0",
                alpha=0.4,
                scale="area",
                trim=FALSE,
                draw_quantiles=c(0.25, 0.5, 0.75)
  )+labs(x=covars[i],
         y="Species"
  )+theme(legend.position="none"
  )+coord_cartesian(xlim=c(xmin,xmax)
  )+theme_minimal(
  )+theme(panel.grid.major=element_line(color="#CCCCCC"),
          panel.grid.minor.x=element_blank())
  
  # ggsave(filename=paste(getwd(),'/EcologicalNichePlots/',covars[i],"_SpeciesPresDensity.png",sep=""),
  #        device="png",
  #        width=600,
  #        height=300,
  #        units="px",
  #        scale=2.5)
  ggsave(filename=paste(getwd(),'/EcologicalNichePlots/',covars[i],"_SpeciesPresDensity_flipped.png",sep=""),
         device="png",
         width=300,
         height=200,
         units="px",
         scale=7)
}




## Compare density curves of covar values experienced across sites ----------------

sites = unique(masterDF$Site)
# sites = c("HAT","GS","BP","BS")
# covars = c("Chl0","EKE0","FSLE0","AEddyDist0","CEddyDist0","SSH0","Temp0","Temp700","Sal0","Sal700",
#           "VelMag0","VelMag700","VelAsp0","VelAsp700")
covars = c("Chl0","FSLE0","SSH0","Temp0","Sal0","VelMag0","VelAsp0")
lineOpts = c("dotdash","solid","dotted","dashed",
             "dotdash","solid","dotted","dashed",
             "dotdash","solid")
masterDF$SiteFac = factor(masterDF$Site,levels=sites)

for (i in 1:length(covars)){
  
  siteInd = which(masterDF$Site%in%sites)
  xmin = min(masterDF[[covars[i]]],na.rm=TRUE)
  xmax = max(masterDF[[covars[i]]],na.rm=TRUE)
  # yrInd = which(masterDF$WeekDate[siteInd]<as.Date("2017-05-01"))
  ggplot(masterDF[siteInd,],aes(x=.data[[covars[i]]])
  )+geom_density(aes(color=SiteFac,linetype=SiteFac)
  )+scale_linetype_manual(values=lineOpts[1:length(sites)]
  )+scale_color_manual(values=c("#f05a43","#db51a4","#de76f5","#9b85f2","#56b4fc","#2fc7cc","#3dd17d","#5aba30","#bbc22f","#d6a10f")
  )+theme(legend.position="none"
  )+coord_cartesian(xlim=c(xmin,xmax)
  )+theme_minimal(
  )+theme(panel.grid.major=element_line(color="#CCCCCC"),
          panel.grid.minor.x=element_blank())
  ggsave(filename=paste(getwd(),'/EcologicalNichePlots/',covars[i],"_Density_1.png",sep=""),device="png",
         width=300,height=150,units="px",scale=7)
  # ggsave(filename=paste(getwd(),'/EcologicalNichePlots/',"ForLegend.png",sep=""),device="png",
  #        width=300,height=300,units="px",scale=7)
  
}

# # sanity check that values are what we expect
# for (i in 1:length(sites)){
#   siteInd = which(masterDF$Site==sites[i])
#   print(paste(sites[i],":",max(masterDF$Chl0[siteInd],na.rm=TRUE)))
#   
# }
## Look at impact of eddies at each site ---------------------
sites = unique(masterDF$Site)
covars = c("Chl0","EKE0","FSLE0","SSH0","Temp0","Sal0",
           "VelMag0","VelAsp0")
ACorrMat = matrix(nrow=length(sites),ncol=length(covars))
CCorrMat = matrix(nrow=length(sites),ncol=length(covars))
ACorrMatAn = matrix(nrow=length(sites),ncol=length(covars))
CCorrMatAn = matrix(nrow=length(sites),ncol=length(covars))
for (i in 1:length(sites)){
  siteInd = which(masterDF$Site==sites[i])
  AplotList = list()
  CplotList = list()
  AplotListAnom = list()
  CplotListAnom = list()
  for (j in 1:length(covars)){
    if (covars[j]!="EKE0"){ # don't calculate anomaly for EKE (it's already an anomaly)
      varMean = rollmean(masterDF[[covars[j]]][siteInd],k=5)
      varAnom = masterDF[[covars[j]]][siteInd] - c(NA,NA,varMean,NA,NA)
    } else {
      varAnom = masterDF[[covars[j]]][siteInd]
    }
    plotDF = data.frame(var=masterDF[[covars[j]]][siteInd],A=masterDF$AEddyDist0[siteInd],C=masterDF$CEddyDist0[siteInd])
    Acorr = cor(plotDF$var,plotDF$A,use="complete.obs")
    ACorrMat[i,j] = round(Acorr,digits=4)
    Ccorr = cor(plotDF$var,plotDF$C,use="complete.obs")
    CCorrMat[i,j] = round(Ccorr,digits=4)
    plotDFAnom = data.frame(anom=varAnom,A=masterDF$AEddyDist0[siteInd],C=masterDF$CEddyDist0[siteInd])
    AcorrAnom = cor(plotDFAnom$anom,plotDFAnom$A,use="complete.obs")
    ACorrMatAn[i,j] = round(AcorrAnom,digits=4)
    CcorrAnom = cor(plotDFAnom$anom,plotDFAnom$C,use="complete.obs")
    CCorrMatAn[i,j] = round(CcorrAnom,digits=4)
    AplotList[[covars[j]]] = ggplot(plotDF
    )+geom_point(aes(x=A,y=var)
    )+labs(y=covars[j],x=NULL
    )+annotate("text",
               label=as.character(round(Acorr,digits=3)),
               x=0.85*max(plotDF$A,na.rm=TRUE),
               y=(0.5*(max(plotDF$var,na.rm=TRUE)-min(plotDF$var,na.rm=TRUE)))+min(plotDF$var,na.rm=TRUE),
               color="blue",
               size=3.5)
    CplotList[[covars[j]]] = ggplot(plotDF
    )+geom_point(aes(x=C,y=var)
    )+labs(y=covars[j],x=NULL
    )+annotate("text",
               label=as.character(round(Ccorr,digits=3)),
               x=0.85*max(plotDF$C,na.rm=TRUE),
               y=(0.5*(max(plotDF$var,na.rm=TRUE)-min(plotDF$var,na.rm=TRUE)))+min(plotDF$var,na.rm=TRUE),
               color="blue",
               size=3.5)
    AplotListAnom[[covars[j]]] = ggplot(plotDFAnom
    )+geom_point(aes(x=A,y=anom)
    )+labs(y=covars[j],x=NULL
    )+annotate("text",
               label=as.character(round(AcorrAnom,digits=3)),
               x=0.85*max(plotDFAnom$A,na.rm=TRUE),
               y=(0.5*(max(plotDFAnom$anom,na.rm=TRUE)-min(plotDFAnom$anom,na.rm=TRUE)))+min(plotDFAnom$anom,na.rm=TRUE),
               color="blue",
               size=3.5)
    CplotListAnom[[covars[j]]] = ggplot(plotDFAnom
    )+geom_point(aes(x=C,y=anom)
    )+labs(y=covars[j],x=NULL
    )+annotate("text",
               label=as.character(round(CcorrAnom,digits=3)),
               x=0.85*max(plotDFAnom$C,na.rm=TRUE),
               y=(0.5*(max(plotDFAnom$anom,na.rm=TRUE)-min(plotDFAnom$anom,na.rm=TRUE)))+min(plotDFAnom$anom,na.rm=TRUE),
               color="blue",
               size=3.5)
  }
  # png(file=paste(getwd(),'/EcologicalNichePlots/',sites[i],'_AEddyInfluence.png',sep=""),height=700,width=600,units="px",res=100)
  # grid.arrange(grobs=AplotList,ncol=2,nrow=5,top=sites[i])
  # while (dev.cur()>1) {dev.off()}
  #
  # png(file=paste(getwd(),'/EcologicalNichePlots/',sites[i],'_CEddyInfluence.png',sep=""),height=700,width=600,units="px",res=100)
  # grid.arrange(grobs=CplotList,ncol=2,nrow=5,top=sites[i])
  # while (dev.cur()>1) {dev.off()}
  #
  # png(file=paste(getwd(),'/EcologicalNichePlots/',sites[i],'_AEddyInfluence_Anom.png',sep=""),height=700,width=600,units="px",res=100)
  # grid.arrange(grobs=AplotListAnom,ncol=2,nrow=5,top=sites[i])
  # while (dev.cur()>1) {dev.off()}
  #
  # png(file=paste(getwd(),'/EcologicalNichePlots/',sites[i],'_CEddyInfluence_Anom.png',sep=""),height=700,width=600,units="px",res=100)
  # grid.arrange(grobs=CplotListAnom,ncol=2,nrow=5,top=sites[i])
  # while (dev.cur()>1) {dev.off()}
}
rownames(ACorrMat) = paste('AEddyDist_',sites,sep="")
rownames(CCorrMat) = paste('CEddyDist_',sites,sep="")
rownames(ACorrMatAn) = paste('AEddyDist_',sites,sep="")
rownames(CCorrMatAn) = paste('CEddyDist_',sites,sep="")
colnames(ACorrMat) = covars
colnames(CCorrMat) = covars
colnames(ACorrMatAn) = covars
colnames(CCorrMatAn) = covars
png(file='EddyDistCorrelations.png',height=600,width=600,units="px",res=75)
par(mfrow=c(2,2))
corrplot(ACorrMat,method='circle')
corrplot(CCorrMat,method='circle')
corrplot(ACorrMatAn,method='circle')
corrplot(CCorrMatAn,method='circle')
while (dev.cur()>1) {dev.off()}

# PCA to compare combinations of conditions targeted by competitive species at co-occupied sites --------------
specs = c("Gm","Gg")
sites = c("HZ","BC","NFC","HAT")
covars = c("Chl0","FSLE0","SSH0","Sal0","Temp0","Sal400","Temp400","AEddyDist0","CEddyDist0")
plotColors = c('#D62A1C','#1C3FD6')
masterDF = data.frame(read.csv('MasterWeeklyDF.csv'))
# remove NAs
badRows = unique(which(is.na(masterDF),arr.ind=TRUE)[,1])
masterDF = masterDF[-badRows,]
masterDF$SiteFac = factor(masterDF$Site,levels=sites)
whichCols= which(colnames(masterDF)%in%covars)
for (i in 1:length(sites)){
  siteInd = which(masterDF$Site==sites[i])
  trainDat = numeric()
  for (k in 1:length(specs)){
    thisSpec = which(masterDF[[specs[k]]][siteInd]>0)
    trainDat = rbind(trainDat,cbind(masterDF[siteInd[thisSpec],whichCols],specs[k]))
  }
  colnames(trainDat)[length(covars)+1] = "Spec"
  pc = prcomp(trainDat[,1:length(whichCols)],center=TRUE,scale=TRUE)
  g <- ggbiplot(pc,
                obs.scale = 1,
                var.scale = 1,
                groups = trainDat$Spec,
                ellipse = TRUE,
                circle = FALSE,
                ellipse.prob = 0.66,
                varname.size=4)
  g = g + scale_color_manual(values=c("#DB51A4","#2FC7CC"))
  g <- g + theme(legend.direction = 'vertical',
                 legend.position = 'right')
  png(filename=paste(getwd(),'/EcologicalNichePlots/',paste(specs,collapse="_"),"_",sites[i],'_PCA.png',sep=""),
      height=600, width=600,units="px",res=100)
  print(g)
  while (dev.cur()>1) {dev.off()}
}