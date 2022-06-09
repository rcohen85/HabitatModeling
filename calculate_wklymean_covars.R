library(ggplot2)
library(zoo)
library(gridExtra)
library(stringr)

data = read.csv("J:/Chpt_3/ModelData/Cuvier_masterDF.csv")
outDir = "J:/Chpt_3/CovarPlots"

covarList = c("Chl0",
              "SSH0",
              "FSLE0",
              "EKE0",
              "AEddyDist0",
              "CEddyDist0",
              "Sal0",
              "Temp0")

sites = c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS")
k=7
startDate = as.Date("2016-05-01",origin="1970-01-01")
endDate = as.Date("2019-04-30",origin="1970-01-01")
dateMarks = c(seq.Date(startDate,endDate,by="year"),endDate)

for (j in 1:length(sites)){
  siteInd = which(!is.na(str_match(data$Site,sites[j])))
    plotDF = data.frame(Date = as.Date(rollmean(data$Date[siteInd],k=k,fill=NA),origin="1970-01-01"))
    for (m in 1:(length(dateMarks)-1)){ # for each year
      plotList = list()
      st = dateMarks[m]
      ed = dateMarks[m+1]
      for (i in 1:length(covarList)){ # plot single layer covars in this year
        plotDF$Covar = rollmean(data[[covarList[i]]][siteInd],k=k,fill=NA)
        plotList[[covarList[i]]] = ggplot(plotDF,aes(x=Date)
        )+geom_line(aes(y=Covar)
        )+scale_x_continuous(breaks=c(seq.Date(from=st,to=ed,by="quarter"))
        )+coord_cartesian(xlim=c(st,ed)
        )+ggtitle(covarList[i]
        )+labs(x=NULL,y=NULL
        )+theme(plot.margin=margin(c(0,43,0,0),unit="pt"))
      }
      png(file=paste(outDir,"/",sites[j],"_",substr(st,1,4),"-",substr(ed,1,4),"_rm.png",sep=""),width=1500,height=1000)
      grid.arrange(grobs=plotList,ncol=1,nrow=length(covarList),top=sites[j])
      while (dev.cur()>1) {dev.off()}
    }
} 
