library(gridExtra)
library(ggplot2)
library(dplyr)
library(rlang)
library(lubridate)
library(stringr)
library(stats)
library(pracma)
library(splines)

inDir = 'J:/Chpt_3/CovarTS'
outDir = 'J:/Chpt_3/CovarPlots'

fileList = dir(inDir,pattern=".Rdata")
# lags = c(7,14,21,28,42,56)

## Clean FSLE data -------------------------
# Aviso FSLE is daily at 1/25the (0.04) degree spatial resolution, re-gridded to 0.08deg
load('J:/Chpt_3/CovarTS/FSLE_TS.Rdata')
FSLE = data.frame(Time=as.Date(masterData.Time[1,],origin="1970-01-01"),HZ0=masterData.Data[1,],
                  OC0=masterData.Data[2,],NC0=masterData.Data[3,],BC0=masterData.Data[4,],
                  WC0=masterData.Data[5,],NFC0=masterData.Data[6,],HAT0=masterData.Data[7,],GS0=masterData.Data[8,],
                  BP0=masterData.Data[9,],BS0=masterData.Data[10,],JAX0=masterData.Data[11,])

sites = colnames(FSLE)

# Plot histograms of data
for (i in 2:12){
  eval(parse(text=(paste(sites[i],' = ggplot(data=FSLE)+geom_histogram(aes(x=',sites[i],'))+labs(x="",title=sites[i])',sep=""))))
}
png(file=paste(outDir,'/',"FSLE_hist.png",sep=""),width = 800, height = 800, units = "px")
grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0,JAX0, ncol=4,nrow=3,top="FSLE")
while (dev.cur()>1) {dev.off()}

# Plot boxplots and calculate quantiles to identify outliers
ggplot(stack(FSLE[,2:12]))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y="FSLE",title="FSLE")
ggsave(paste(outDir,'/',"FSLE_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)

q25 = quantile(stack(FSLE[,2:12])$values,probs=0.25)
q75 = quantile(stack(FSLE[,2:12])$values,probs=0.75)
iqr = q75-q25

# Remove outliers (make NA, interpolate later)
# for (i in 1:11){
# which_outliers = which(FSLE[,i]<(q25-(1.5*iqr)) | FSLE[,i]>0)
# FSLE[which_outliers,i] = NA}

# Check for missing dates
timeDiff = diff(as.numeric(FSLE$Time))
any(timeDiff>1)

# Interpolate missing values
for (i in 2:12){
  skippedBins = which(is.na(FSLE[,i])) # missing data
  missFSLE = apply(cbind(FSLE[skippedBins-1,i],FSLE[skippedBins+1,i]),MARGIN=1,mean)
  FSLE[skippedBins,i] = missFSLE
}

# # Create time lagged vectors
# startInd = which(FSLE$Time==as.Date('2016-05-01',origin="1970-01-01"))
# fullLength = length(FSLE$Time)
# for (i in 2:12){
#   for (k in lags){
#     lagInd = startInd-k
#     eval(parse(text=paste('FSLE$',sites[i],'Lag',as.character(k),' = NA',sep="")))
#     eval(parse(text=paste('FSLE$',sites[i],'Lag',as.character(k),'[startInd:',fullLength,'] = FSLE$',sites[i],'[lagInd:(',fullLength,'-',as.character(k),')]',sep="")))
#   }
# }
write.csv(FSLE,file=paste(inDir,'/','FSLE_TS.csv',sep=""),row.names=FALSE)

## Clean HYCOM data --------------------------------
# Downloaded HYCOM data are daily at 2/25th (0.08) degrees irregular spatial resolution, re-gridded to 0.08deg
vars = c('SSH','Salinity','Temperature','VelocityMag','VelocityAsp','EKE')
lon = "ES"

for (j in 1:length(vars)){
  
  if (vars[j]=="SSH"){
    
    TSind = which(!is.na(str_match(fileList,vars[j])) & !is.na(str_match(fileList,lon)))
    load(paste(inDir,'/',fileList[TSind],sep=""))
    sites = c('HZ0','OC0','NC0','BC0','WC0','NFC0','HAT0','GS0','BP0','BS0','JAX0')
    
    SSHDF = data.frame(t(masterData.Data))
    colnames(SSHDF) = sites
    
    # Histograms
    for (i in 1:11){
      eval(parse(text=(paste(sites[i],' = ggplot(data=SSHDF)+geom_histogram(aes(x=',sites[i],'))+labs(x="m",y="",title=sites[i])',sep=""))))
    }
    png(file=paste(outDir,'/',"SSH_hist.png",sep=""),width = 800, height = 800, units = "px")
    grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0,JAX0, ncol=4,nrow=3,top="SSH")
    while (dev.cur()>1) {dev.off()}
    
    # Plot boxplots and calculate quantiles to identify outliers
    ggplot(stack(SSHDF))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y="m",title="SSH")
    ggsave(paste(outDir,'/',"SSH_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)
    
    q25 = quantile(stack(SSHDF)$values,probs=0.25)
    q75 = quantile(stack(SSHDF)$values,probs=0.75)
    iqr = q75-q25
    
    # Remove outliers?
    for (i in 1:11){
      q = round(SSHDF[,i],digits=10)
      SSHDF[q==0,i] = NA
    }
    
    # Interpolate missing dates and NAs
    fullSSH = data.frame(Time=seq.Date(from=as.Date("2016-02-01",origin="1970-01-01"),to=as.Date("2019-04-30",origin="1970-01-01"),by=1))
    for (i in 1:11){
      datBins = which(!is.na(SSHDF[,i]))
      eval(parse(text=paste('fullSSH$',sites[i],'=NA',sep="")))
      eval(parse(text=paste('fullSSH$',sites[i],' = (approx(x=masterData.Time[datBins],y=SSHDF[datBins,i],xout=fullSSH[,1],method="linear"))$y',sep="")))
    }
    
    # # Create time lagged vectors
    # startInd = which(fullSSH$Time==as.Date('2016-05-01',origin="1970-01-01"))
    # fullLength = length(fullSSH$Time)
    # for (i in 1:10){
    #   for (k in lags){
    #     lagInd = startInd-k
    #     eval(parse(text=paste('fullSSH$',sites[i],'Lag',as.character(k),' = NA',sep="")))
    #     eval(parse(text=paste('fullSSH$',sites[i],'Lag',as.character(k),'[startInd:',fullLength,'] = fullSSH$',sites[i],'[lagInd:(',fullLength,'-',as.character(k),')]',sep="")))
    #   }
    # }
    write.csv(fullSSH,file=paste(inDir,'/','SSH_TS.csv',sep=""),row.names=FALSE)
    
  } else { # Multi-depth Variables
    
    TSind = which(!is.na(str_match(fileList,paste(vars[j],'_Profiles',sep=""))) & !is.na(str_match(fileList,lon)))
    load(paste(inDir,'/',fileList[TSind],sep=""))
    outDFs = c('fullSalinity','fullTemperature','fullVelocityMag','fullVelocityAsp','fullEKE')
    depths = c(0,100,200,300,400,500,600,700,800) # desired depth layers
    units = c("PPT",paste(intToUtf8(176),"C",sep=""),"m/s",intToUtf8(176),paste(expression('cm'^2),'/',expression('s'^2),sep=""))
    
    eval(parse(text=paste('full',vars[j],' = data.frame(Time=seq.Date(from=as.Date("2016-02-01",origin="1970-01-01"),to=as.Date("2019-04-30",origin="1970-01-01"),by=1))',sep="")))
    
    for (l in 1:length(depths)){
      
      eval(parse(text=paste('Temp = data.frame(HZ',depths[l],'=HZProfile[l,],
                            OC',depths[l],'=OCProfile[l,],
                            NC',depths[l],'=NCProfile[l,],
                            BC',depths[l],'=BCProfile[l,],
                            WC',depths[l],'=WCProfile[l,],
                            NFC',depths[l],'=NFCProfile[l,],
                            HAT',depths[l],'=HATProfile[l,],
                            GS',depths[l],'=GSProfile[l,],
                            BP',depths[l],'=BPProfile[l,],
                            BS',depths[l],'=BSProfile[l,],
                            JAX',depths[l],'=JAXProfile[l,])',sep="")))
      sites = colnames(Temp)
      siteName = str_replace(sites,as.character(depths[l]),"")
      
      for (i in 1:11){
        eval(parse(text=paste('Temp$',sites[i],'[is.nan(Temp$',sites[i], ')] = NA',sep="")))
        eval(parse(text=paste('Temp$',sites[i],'[Temp$',sites[i],'==0] = NA',sep="")))
      }
      
      # Histograms
      for (i in 1:11){
        eval(parse(text=(paste(siteName[i],' = ggplot(data=Temp)+geom_histogram(aes(x=',
                               sites[i],'))+labs(x=units[j-1],title=sites[i])',sep=""))))
      }
      png(file=paste(outDir,'/',vars[j],'_',as.character(depths[l]),"_hist_ES.png",sep=""),width = 800, height = 800, units = "px")
      grid.arrange(HZ,NC,OC,BC,WC,NFC,HAT,GS,BP,BS,JAX, ncol=4,nrow=3,top=paste(vars[j],' at ',as.character(depths[l]),'m',sep=""))
      while (dev.cur()>1) {dev.off()}
      
      # Plot boxplots and calculate quantiles to identify outliers
      ggplot(stack(Temp))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y=units[j-1],title=paste(vars[j],' at ',as.character(depths[l]),'m',sep=""))
      ggsave(paste(outDir,'/',vars[j],'_',as.character(depths[l]),"_boxplot_ES.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)
      
      # Remove outliers?
      q25 = quantile(stack(Temp)$values,probs=0.25,na.rm=TRUE)
      q75 = quantile(stack(Temp)$values,probs=0.75,na.rm=TRUE)
      iqr = q75-q25
      upperThresh = q75+(1.5*iqr)
      lowerThresh = q25-(1.5*iqr)
      
      # Interpolate missing dates and NAs (can't interpolate JAX at most depths, too little data!)
      if (l==1){    
        for (i in 1:11){
          datBins = which(!is.na(Temp[,i]))
          eval(parse(text=paste('full',vars[j],'$',sites[i],'=NA',sep="")))
          eval(parse(text=paste('full',vars[j],'$',sites[i],' = (approx(x=masterData.Time[datBins],y=Temp[datBins,i],xout=full',vars[j],'[,1],method="linear"))$y',sep="")))
        }
      } else {
        for (i in 1:10){
          datBins = which(!is.na(Temp[,i]))
          eval(parse(text=paste('full',vars[j],'$',sites[i],'=NA',sep="")))
          eval(parse(text=paste('full',vars[j],'$',sites[i],' = (approx(x=masterData.Time[datBins],y=Temp[datBins,i],xout=full',vars[j],'[,1],method="linear"))$y',sep="")))
        }
        # add (likely hole-y) JAX data without interpolation
        eval(parse(text=paste('full',vars[j],'$',sites[i+1],'=NA',sep="")))
        eval(parse(text=paste('putWhere = match(masterData.Time,full',vars[j],'$Time)',sep="")))
        eval(parse(text=paste('full',vars[j],'$',sites[i+1],'[putWhere[!is.na(putWhere)]]=Temp[-which(is.na(putWhere)),i+1]',sep="")))
      }
      
    }
    
    if (eval(parse(text=paste('any(full',vars[j],'<0,na.rm=TRUE)',sep="")))){
      return('WARNING: SPURIOUS NEGATIVE VALUES CHECK DATA')
    }
    
    # # Create time lagged vectors
    # eval(parse(text=paste('startInd = which(full',vars[j],'$Time==as.Date("2016-05-01",origin="1970-01-01"))',sep="")))
    # eval(parse(text=paste('fullLength = length(full',vars[j],'$Time)',sep="")))
    # eval(parse(text=paste('fullColNames = colnames(full',vars[j],')',sep="")))
    # for (i in 2:length(fullColNames)){
    #   for (k in lags){
    #     lagInd = startInd-k
    #     eval(parse(text=paste('full',vars[j],'$',fullColNames[i],'Lag',as.character(k),' = NA',sep="")))
    #     eval(parse(text=paste('full',vars[j],'$',fullColNames[i],'Lag',as.character(k),'[startInd:',fullLength,'] = full',vars[j],'$',fullColNames[i],'[lagInd:(',fullLength,'-',as.character(k),')]',sep="")))
    #   }
    # }
    
    # Save data as .csv
    eval(parse(text=paste("write.csv(full",vars[j],",file=\"J:/Chpt_3/CovarTS/",vars[j],"_TS.csv\",row.names=FALSE)",sep="")))
    
  }
}

## Clean Chla Data  --------------------------------------
# Downloaded Chla data are daily at 1/24th (0.04166667) degree spatial resolution, effectively re-gridded to 0.08deg
load('J:/Chpt_3/CovarTS/Chl_TS_4Averaged.Rdata')
Chl = data.frame(Time=as.Date(masterData.Time,origin="1970-01-01"),HZ0=masterData.Data[1,],
                 OC0=masterData.Data[2,],NC0=masterData.Data[3,],BC0=masterData.Data[4,],
                 WC0=masterData.Data[5,],NFC0=masterData.Data[6,],HAT0=masterData.Data[7,],GS0=masterData.Data[8,],
                 BP0=masterData.Data[9,],BS0=masterData.Data[10,],JAX0=masterData.Data[11,])

sites = colnames(Chl)

# Plot histograms
for (i in 2:12){
  eval(parse(text=(paste(sites[i],' = ggplot(data=Chl)+geom_histogram(aes(x=',sites[i],'))+labs(x=expression("mg/m"^3),title=sites[i])',sep=""))))
}
png(file=paste(outDir,'/',"Chl_hist.png",sep=""),width = 800, height = 800, units = "px")
grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0,JAX0, ncol=4,nrow=3,top="Chl")
while (dev.cur()>1) {dev.off()}

# Plot boxplots and calculate quantiles to identify outliers
ggplot(stack(Chl[,2:12]))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y=expression("mg/m"^3),title="Chl")
ggsave(paste(outDir,'/',"Chl_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)

q25 = quantile(stack(Chl)$values,probs=0.25,na.rm=TRUE)
q75 = quantile(stack(Chl)$values,probs=0.75,na.rm=TRUE)
iqr = q75-q25

# Remove outliers?

# Check for missing dates
timeDiff = diff(as.numeric(Chl$Time))
any(timeDiff>1)

# Interpolate missing values (many missing, need to interpolate whole TS at each
# site, instead of just interpolating across occasional gaps)
for (i in 2:12){
  datBins = which(!is.na(Chl[,i]))
  Chl[,i] = (approx(x=Chl[datBins,1],y=Chl[datBins,i],xout=Chl[,1],method="linear"))$y
}

# # Create time lagged vectors
# startInd = which(Chl$Time==as.Date('2016-05-01',origin="1970-01-01"))
# fullLength = length(Chl$Time)
# for (i in 2:12){
#   for (k in lags){
#     lagInd = startInd-k
#     eval(parse(text=paste('Chl$',sites[i],'Lag',as.character(k),' = NA',sep="")))
#     eval(parse(text=paste('Chl$',sites[i],'Lag',as.character(k),'[startInd:',fullLength,'] = Chl$',sites[i],'[lagInd:(',fullLength,'-',as.character(k),')]',sep="")))
#   }
# }
write.csv(Chl,file=paste(inDir,'/','Chl_TS.csv',sep=""),row.names=FALSE)

## Clean eddy data? ---------------------------
load('J:/Chpt_3/CovarTS/EddyDist_TS.Rdata')
eddyA = data.frame(Time=as.Date(masterData.Time[1,],origin="1970-01-01"),HZ0=masterData.eddyDistA[1,],
                 OC0=masterData.eddyDistA[2,],NC0=masterData.eddyDistA[3,],BC0=masterData.eddyDistA[4,],
                 WC0=masterData.eddyDistA[5,],NFC0=masterData.eddyDistA[6,],HAT0=masterData.eddyDistA[7,],
                 GS0=masterData.eddyDistA[8,],BP0=masterData.eddyDistA[9,],BS0=masterData.eddyDistA[10,],
                 JAX0=masterData.eddyDistA[11,])
eddyC = data.frame(Time=as.Date(masterData.Time[1,],origin="1970-01-01"),HZ0=masterData.eddyDistC[1,],
                   OC0=masterData.eddyDistC[2,],NC0=masterData.eddyDistC[3,],BC0=masterData.eddyDistC[4,],
                   WC0=masterData.eddyDistC[5,],NFC0=masterData.eddyDistC[6,],HAT0=masterData.eddyDistC[7,],
                   GS0=masterData.eddyDistC[8,],BP0=masterData.eddyDistC[9,],BS0=masterData.eddyDistC[10,],
                   JAX0=masterData.eddyDistC[11,])

sites = colnames(eddyA)

# Plot histograms
for (i in 2:12){
  eval(parse(text=(paste(sites[i],' = ggplot(data=eddyA)+geom_histogram(aes(x=',sites[i],'))+labs(x=expression("km"),title=sites[i])',sep=""))))
}
png(file=paste(outDir,'/',"AEddyDist_hist.png",sep=""),width = 800, height = 800, units = "px")
grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0,JAX0, ncol=4,nrow=3,top="Distance to Anticyclonic Eddies")
while (dev.cur()>1) {dev.off()}

for (i in 2:12){
  eval(parse(text=(paste(sites[i],' = ggplot(data=eddyC)+geom_histogram(aes(x=',sites[i],'))+labs(x=expression("km"),title=sites[i])',sep=""))))
}
png(file=paste(outDir,'/',"CEddyDist_hist.png",sep=""),width = 800, height = 800, units = "px")
grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0,JAX0, ncol=4,nrow=3,top="Distance to Cyclonic Eddies")
while (dev.cur()>1) {dev.off()}

# Plot boxplots and calculate quantiles to identify outliers
ggplot(stack(eddyA[,2:12]))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y=expression("km"),title="Distance to Anticyclonic Eddies")
ggsave(paste(outDir,'/',"AEddyDist_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)

ggplot(stack(eddyC[,2:12]))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y=expression("km"),title="Distance to Cyclonic Eddies")
ggsave(paste(outDir,'/',"CEddyDist_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)

# Check for missing dates
timeDiff = diff(as.numeric(eddyA$Time))
any(timeDiff>1)

# doesn't make sense to interpolate this time series

# save
write.csv(eddyA,file=paste(inDir,'/','AEddyDist_TS.csv',sep=""),row.names=FALSE)
write.csv(eddyC,file=paste(inDir,'/','CEddyDist_TS.csv',sep=""),row.names=FALSE)
