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

fileList = dir(inDir)

## Clean FSLE data -------------------------
# Aviso FSLE is daily at 1/25the (0.04) degree spatial resolution, regridded to 0.08deg
# load('J:/Chpt_3/CovarTS/FSLE_TS.Rdata')
# FSLE = data.frame(Time=as.Date(masterData.Time[1,],origin="1970-01-01"),HZ0=masterData.Data[1,],
#                   OC0=masterData.Data[2,],NC0=masterData.Data[3,],BC0=masterData.Data[4,],
#                      WC0=masterData.Data[5,],NFC0=masterData.Data[6,],HAT0=masterData.Data[7,],GS0=masterData.Data[8,],
#                      BP0=masterData.Data[9,],BS0=masterData.Data[10,],JAX0=masterData.Data[11,])
# 
# sites = colnames(FSLE)
# 
# # Plot histograms of data
# for (i in 2:12){
#   eval(parse(text=(paste(sites[i],' = ggplot(data=FSLE)+geom_histogram(aes(x=',sites[i],'))+labs(x="",title=sites[i])',sep=""))))
# }
# png(file=paste(outDir,'/',"FSLE_hist.png",sep=""),width = 800, height = 800, units = "px")
# grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0,JAX0, ncol=4,nrow=3,top="FSLE")
# while (dev.cur()>1) {dev.off()}
# 
# # Plot boxplots and calculate quantiles to identify outliers
# ggplot(stack(FSLE[,2:12]))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y="FSLE",title="FSLE")
# ggsave(paste(outDir,'/',"FSLE_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)
# 
# q25 = quantile(stack(FSLE[,2:12])$values,probs=0.25)
# q75 = quantile(stack(FSLE[,2:12])$values,probs=0.75)
# iqr = q75-q25
# 
# # Remove outliers (make NA, interpolate later)
# # for (i in 1:11){
# # which_outliers = which(FSLE[,i]<(q25-(1.5*iqr)) | FSLE[,i]>0)
# # FSLE[which_outliers,i] = NA}
# 
# # Check for missing dates
# timeDiff = diff(as.numeric(FSLE$Time))
# any(timeDiff>1)
# 
# # Interpolate missing values
# for (i in 2:12){
#   skippedBins = which(is.na(FSLE[,i])) # missing data
#   missFSLE = apply(cbind(FSLE[skippedBins-1,i],FSLE[skippedBins+1,i]),MARGIN=1,mean)
#   FSLE[skippedBins,i] = missFSLE
# }
# 
# # Create time lagged vectors
# startInd = which(FSLE$Time==as.Date('2016-05-01',origin="1970-01-01"))
# fullLength = length(FSLE$Time)
# for (i in 2:12){
#   for (k in 1:3){
#     if (k==1){
#       lag = 7
#       lagInd = startInd-lag
#     } else if (k==2){
#       lag=14
#       lagInd = startInd-lag
#     } else if (k==3){
#       lag=21
#       lagInd = startInd-lag
#     }
#   eval(parse(text=paste('FSLE$',sites[i],'Lag',lag,' = NA',sep="")))
#   eval(parse(text=paste('FSLE$',sites[i],'Lag',lag,'[startInd:',fullLength,'] = FSLE$',sites[i],'[lagInd:(',fullLength,'-lag)]',sep="")))
#   }
# }
# write.csv(FSLE,file=paste(outDir,'/','FSLE_TS.csv',sep=""),row.names=FALSE)

## Clean HYCOM data --------------------------------
# Downloaded HYCOM data are daily at 2/25th (0.08) degrees spatial resolution
vars = c('SSH','Salinity','Temperature','VelocityMag','VelocityAsp')
for (j in 1:length(vars)){

  if (vars[j]=="SSH"){

    TSind = which(!is.na(str_match(fileList,vars[j])))
    load(paste(inDir,'/',fileList[TSind],sep=""))
    sites = c('HZ0','OC0','NC0','BC0','WC0','NFC0','HAT0','GS0','BP0','BS0')

    SSHDF = data.frame(t(masterData.Data))
    colnames(SSHDF) = sites

    # Histograms
    for (i in 1:10){
      eval(parse(text=(paste(sites[i],' = ggplot(data=SSHDF)+geom_histogram(aes(x=',sites[i],'))+labs(x="m",y="",title=sites[i])',sep=""))))
    }
    png(file=paste(outDir,'/',"SSH_hist.png",sep=""),width = 800, height = 800, units = "px")
    grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0, ncol=4,nrow=3,top="SSH")
    while (dev.cur()>1) {dev.off()}

    # Plot boxplots and calculate quantiles to identify outliers
    ggplot(stack(SSHDF))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y="m",title="SSH")
    ggsave(paste(outDir,'/',"SSH_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)

    q25 = quantile(stack(SSHDF)$values,probs=0.25)
    q75 = quantile(stack(SSHDF)$values,probs=0.75)
    iqr = q75-q25

    # Remove outliers?
    # for (i in 1:11){
    #   q = round(plotDF[,i],digits=10)
    #   plotDF[q==0,i] = NA
    # }

    # Check for missing dates
    timeDiff = diff(masterData.Time)
    skipped = which(timeDiff>1)

    # Interpolate missing dates
    fullSSH = data.frame(Time=seq.Date(from=as.Date("2016-02-01",origin="1970-01-01"),to=as.Date("2019-04-30",origin="1970-01-01"),by=1))
    for (i in 1:10){
      # fixedHYCOMTime = floor_date(as_date((HYCOM$Time[[i]][[1]][1:chopInd]-719529),origin="1970-01-01"),unit="day")
      putWhere1 = match(masterData.Time,fullSSH$Time)
      eval(parse(text=paste('fullSSH$',sites[i],'=NA',sep="")))
      eval(parse(text=paste('fullSSH$',sites[i],'[putWhere1[!is.na(putWhere1)]]=SSHDF$',sites[i],'[-which(is.na(putWhere1))]',sep="")))
      missSSH = apply(cbind(SSHDF[skipped,i],SSHDF[skipped+1,i]),MARGIN=1,mean)
      missTimes = masterData.Time[skipped]+1
      putWhere2 = match(skipped,fullSSH$Time) # indices where data is missing
      missSSH = missSSH[-which(is.na(putWhere2))]
      putWhere2 = putWhere2[-which(is.na(putWhere2))]
      eval(parse(text=paste('fullSSH$',sites[i],'[putWhere2] = missSSH',sep="")))
    }

    # Check for NAs and interpolate
    for (i in 1:10){
      skippedBins = which(is.na(fullSSH[,i])) # missing data
      missSSH = apply(cbind(fullSSH[skippedBins-1,i],fullSSH[skippedBins+1,i]),MARGIN=1,mean)
      fullSSH[skippedBins,i] = missSSH
    }

    # Create time lagged vectors
    startInd = which(fullSSH$Time==as.Date('2016-05-01',origin="1970-01-01"))
    fullLength = length(fullSSH$Time)
    for (i in 1:10){
      for (k in 1:3){
        if (k==1){
          lag = 7
          lagInd = startInd-lag
        } else if (k==2){
          lag=14
          lagInd = startInd-lag
        } else if (k==3){
          lag=21
          lagInd = startInd-lag
        }
        eval(parse(text=paste('fullSSH$',sites[i],'Lag',lag,' = NA',sep="")))
        eval(parse(text=paste('fullSSH$',sites[i],'Lag',lag,'[startInd:',fullLength,'] = fullSSH$',sites[i],'[lagInd:(',fullLength,'-lag)]',sep="")))
      }
    }
    write.csv(fullSSH,file=paste(outDir,'/','SSH_TS.csv',sep=""),row.names=FALSE)
  } else { # Multi-depth Variables

    TSind = which(!is.na(str_match(fileList,paste(vars[j],'_Profiles',sep=""))))
    load(paste(inDir,'/',fileList[TSind],sep=""))
    outDFs = c('fullSalinity','fullTemperature','fullVelocityMag','fullVelocityAsp')
    depths = c(0,100,200,300,400,500,600,700,800) # desired depth layers
    units = c("PPT",paste(intToUtf8(176),"C",sep=""),"m/s",'m/s')

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
                            BS',depths[l],'=BSProfile[l,])',sep="")))
      sites = colnames(Temp)
      siteName = str_replace(sites,as.character(depths[l]),"")

      for (i in 1:10){
        eval(parse(text=paste('Temp$',sites[i],'[is.nan(Temp$',sites[i], ')] = NA',sep="")))
        eval(parse(text=paste('Temp$',sites[i],'[Temp$',sites[i],'==0] = NA',sep="")))
      }

      # Histograms
      for (i in 1:10){
        eval(parse(text=(paste(siteName[i],' = ggplot(data=Temp)+geom_histogram(aes(x=',
                               sites[i],'))+labs(x=units[j],title=sites[i])',sep=""))))
      }
      png(file=paste(outDir,'/',vars[j],'_',as.character(depths[l]),"_hist.png",sep=""),width = 800, height = 800, units = "px")
      grid.arrange(HZ,NC,OC,BC,WC,NFC,HAT,GS,BP,BS, ncol=4,nrow=3,top=paste(vars[j],' at ',as.character(depths[l]),'m',sep=""))
      while (dev.cur()>1) {dev.off()}

      # Plot boxplots and calculate quantiles to identify outliers
      ggplot(stack(Temp))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y="m",title=paste(vars[j]),' at ',as.character(depths[l]),'m',sep="")
      ggsave(paste(outDir,'/',vars[j],'_',as.character(depths[l]),"_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)

      # Remove outliers?
      q25 = quantile(stack(Temp)$values,probs=0.25,na.rm=TRUE)
      q75 = quantile(stack(Temp)$values,probs=0.75,na.rm=TRUE)
      iqr = q75-q25
      upperThresh = q75+(1.5*iqr)
      lowerThresh = q25-(1.5*iqr)

      # Check for missing dates
        timeDiff = diff(masterData.Time)
        skipped = which(timeDiff>1)

      # Interpolate missing dates
        if (!isempty(skipped)){
          for (i in 1:10){
            # fixedHYCOMTime = floor_date(as_date((HYCOM$Time[[i]][[1]][1:chopInd]-719529),origin="1970-01-01"),unit="day")
            eval(parse(text=paste('putWhere1 = match(masterData.Time,full',vars[j],'$Time)',sep="")))
            eval(parse(text=paste('full',vars[j],'$',sites[i],'=NA',sep="")))
            eval(parse(text=paste('full',vars[j],'$',sites[i],'[putWhere1[!is.na(putWhere1)]]=Temp$',sites[i],'[-which(is.na(putWhere1))]',sep="")))
            missData = apply(cbind(Temp[skipped,i],Temp[skipped+1,i]),MARGIN=1,mean)
            missTimes = masterData.Time[skipped]+1
            eval(parse(text=paste('putWhere2 = match(missTimes,full',vars[j],'$Time)',sep=""))) # indices where data is missing
            missData = missData[-which(is.na(putWhere2))]
            putWhere2 = putWhere2[-which(is.na(putWhere2))]
            eval(parse(text=paste('full',vars[j],'$',sites[i],'[putWhere2] = missData',sep="")))

          }
        } else {
          for (i in 1:10){
            eval(parse(text=paste('putWhere1 = match(masterData.Time,full',vars[j],'$Time)',sep="")))
            eval(parse(text=paste('full',vars[j],'$',sites[i],'=NA',sep="")))
            eval(parse(text=paste('full',vars[j],'$',sites[i],'[putWhere1[!is.na(putWhere1)]]=Temp$',sites[i],'[-which(is.na(putWhere1))]',sep="")))
          }
        }

    }


    # Check for NAs and linearly interpolate across gaps
    for (i in 2:91){
      eval(parse(text=paste('skippedBins = which(is.na(full',vars[j],'[,i]))',sep=""))) # missing data
      if(any(diff(skippedBins)==1)){
        eval(parse(text=paste('datBins = setdiff(1:dim(full',vars[j],')[1],skippedBins)',sep="")))
        eval(parse(text=paste('missDatNS = approx(x=full',vars[j],'[datBins,1],y=full',vars[j],'[datBins,i],xout=full',vars[j],'[,1],method=\"linear\")',sep="")))
      } else {
      eval(parse(text=paste('missData = apply(cbind(full',vars[j],'[skippedBins-1,i],full',vars[j],'[skippedBins+1,i]),MARGIN=1,mean)',sep="")))
      eval(parse(text=paste('full',vars[j],'[skippedBins,i] = missData',sep="")))
      }
    }

    if (eval(parse(text=paste('any(full',vars[j],'<0)',sep="")))){
      return('WARNING: SPURIOUS NEGATIVE VALUES CHECK DATA')
    }


    # Create time lagged vectors
    eval(parse(text=paste('startInd = which(full',vars[j],'$Time==as.Date("2016-05-01",origin="1970-01-01"))',sep="")))
    eval(parse(text=paste('fullLength = length(full',vars[j],'$Time)',sep="")))
    eval(parse(text=paste('fullColNames = colnames(full',vars[j],')',sep="")))
    for (i in 2:length(fullColNames)){
      for (k in 1:3){
        if (k==1){
          lag = 7
          lagInd = startInd-lag
        } else if (k==2){
          lag=14
          lagInd = startInd-lag
        } else if (k==3){
          lag=21
          lagInd = startInd-lag
        }
        # baseData = paste(sites[i],'0',sep="")
        # thisDepth =
        eval(parse(text=paste('full',vars[j],'$',fullColNames[i],'Lag',lag,' = NA',sep="")))
        eval(parse(text=paste('full',vars[j],'$',fullColNames[i],'Lag',lag,'[startInd:',fullLength,'] = full',vars[j],'$',fullColNames[i],'[lagInd:(',fullLength,'-lag)]',sep="")))
      }
    }

    # Save data as .csv
    eval(parse(text=paste("write.csv(full",vars[j],",file=\"J:/Chpt_3/CovarPlots/",vars[j],"_TS.csv\",row.names=FALSE)",sep="")))

  }

}
## Clean Chla Data  --------------------------------------
# Downloaded Chla data are daily at 1/24th (0.04166667) degree spatial resolution
# load('J:/Chpt_3/CovarTS/Chl_TS_4Averaged.Rdata')
# Chl = data.frame(Time=as.Date(masterData.Time,origin="1970-01-01"),HZ0=masterData.Data[1,],
#                  OC0=masterData.Data[2,],NC0=masterData.Data[3,],BC0=masterData.Data[4,],
#                   WC0=masterData.Data[5,],NFC0=masterData.Data[6,],HAT0=masterData.Data[7,],GS0=masterData.Data[8,],
#                   BP0=masterData.Data[9,],BS0=masterData.Data[10,],JAX0=masterData.Data[11,])
# 
# sites = colnames(Chl)
# 
# # Plot histograms
# for (i in 2:12){
#   eval(parse(text=(paste(sites[i],' = ggplot(data=Chl)+geom_histogram(aes(x=',sites[i],'))+labs(x=expression("mg/m"^3),title=sites[i])',sep=""))))
# }
# png(file=paste(outDir,'/',"Chl_hist.png",sep=""),width = 800, height = 800, units = "px")
# grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0,JAX0, ncol=4,nrow=3,top="Chl")
# while (dev.cur()>1) {dev.off()}
# 
# # Plot boxplots and calculate quantiles to identify outliers
# ggplot(stack(Chl[,2:12]))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y=expression("mg/m"^3),title="Chl")
# ggsave(paste(outDir,'/',"Chl_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)
# 
# q25 = quantile(stack(Chl)$values,probs=0.25,na.rm=TRUE)
# q75 = quantile(stack(Chl)$values,probs=0.75,na.rm=TRUE)
# iqr = q75-q25
# 
# # Remove outliers?
# 
# # Check for missing dates
# timeDiff = diff(as.numeric(Chl$Time))
# any(timeDiff>1)
# 
# # Interpolate missing values (many missing, need to interpolate whole TS at each
# # site, instead of just interpolating across occasional gaps)
# for (i in 2:12){
#   datBins = which(!is.na(Chl[,i]))
#   
#   missChl = approx(x=Chl[datBins,1],y=Chl[datBins,i],xout=Chl[,1],method="linear")
# 
#   Chl[,i] = missChlNS
# }
# 
# # Create time lagged vectors
# startInd = which(Chl$Time==as.Date('2016-05-01',origin="1970-01-01"))
# fullLength = length(Chl$Time)
# for (i in 2:12){
#   for (k in 1:3){
#     if (k==1){
#       lag = 7
#       lagInd = startInd-lag
#     } else if (k==2){
#       lag=14
#       lagInd = startInd-lag
#     } else if (k==3){
#       lag=21
#       lagInd = startInd-lag
#     }
#     eval(parse(text=paste('Chl$',sites[i],'Lag',lag,' = NA',sep="")))
#     eval(parse(text=paste('Chl$',sites[i],'Lag',lag,'[startInd:',fullLength,'] = Chl$',sites[i],'[lagInd:(',fullLength,'-lag)]',sep="")))
#   }
# }
# write.csv(Chl,file=paste(outDir,'/','Chl_TS.csv',sep=""),row.names=FALSE)

## Make modeling csv files ---------------------------

