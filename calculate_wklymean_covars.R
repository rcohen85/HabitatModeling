library(ggplot2)
library(zoo)
library(gridExtra)

data = read.csv("E:/ModelingCovarData/Master_DFs/Cuvier_masterDF.csv")
outDir = "E:/ModelingCovarData/CovarTS_Plots"

covarList = c("Chl0",
              "SSH0",
              "FSLE0",
              "EKE0",
              "AEddyDist",
              "CEddyDist",
              "Sal0",
              "Sal700",
              "Temp0",
              "Temp700",
              "VelMag0",
              "VelMag700",
              "VelAsp0",
              "VelAsp700",
              "GSLat",
              "GSDist")

sites = c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS")
k=7
startDate = as.Date("2016-05-01",origin="1970-01-01")
endDate = as.Date("2019-04-30",origin="1970-01-01")
dateMarks = c(seq.Date(startDate,endDate,by="year"),endDate)

for (j in 1:length(sites)){
  siteInd = which(!is.na(str_match(data$Site,sites[j])))
  # if (sum(data$Pres[siteInd]>0)>25){
    plotDF = data.frame(Pres = rollmean(data$Pres[siteInd],k=k,fill=NA),
                        Date = as.Date(rollmean(data$Date[siteInd],k=k,fill=NA),origin="1970-01-01"))
    for (m in 1:(length(dateMarks)-1)){ # for each year
      # plotList = list()
      st = dateMarks[m]
      ed = dateMarks[m+1]
      # plotList[['Pres']] = ggplot(plotDF,aes(x=Date) # plot presence in this year
      # )+geom_line(aes(y=Pres)
      # )+scale_x_continuous(breaks=c(seq.Date(from=st,to=ed,by="quarter"))
      # )+coord_cartesian(xlim=c(st,ed)
      # )+ggtitle('Presence'
      # )+labs(x=NULL,y=NULL
      # )+theme(plot.margin=margin(c(0,43,0,0),unit="pt"))
      for (i in 1:length(covarList)){ # plot single layer covars in this year
        plotDF$Covar = rollmean(data[[covarList[i]]][siteInd],k=k,fill=NA)
        
        png(filename = paste(outDir,"/",covarList[i],"_",sites[j],"_",substr(st,1,4),"-",substr(ed,1,4),"_rm.png",sep=""))
        ggplot(plotDF,aes(x=Date)
        )+geom_line(aes(y=Covar)
        )+scale_x_continuous(breaks=c(seq.Date(from=st,to=ed,by="quarter"))
        )+coord_cartesian(xlim=c(st,ed)
        )+ggtitle(covarList[i]
        )+labs(x=NULL,y=NULL
        )+theme(plot.margin=margin(c(0,43,0,0),unit="pt"))
        while (dev.cur()>1) {dev.off()}
      }
      # grid.arrange(grobs=plotList,ncol=1,nrow=length(covarList)+1,top=paste('Cuvier at',sites[j]))
    }
  # }
} 
