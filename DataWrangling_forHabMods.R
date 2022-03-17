library(gridExtra)
library(ggplot2)
library(dplyr)
library(rlang)
library(lubridate)
library(stringr)
library(stats)
library(splines2)

outDir = 'J:/Chpt_3/CovarPlots'

## Clean FSLE data -------------------------
# load('J:/Chpt_3/FSLE_TS.Rdata')
# FSLE = data.frame(Time=as.Date(masterData.Time[1,],origin="1970-01-01"),HZ0=masterData.Fsle[1,],OC0=masterData.Fsle[2,],NC0=masterData.Fsle[3,],BC0=masterData.Fsle[4,],
#                      WC0=masterData.Fsle[5,],NFC0=masterData.Fsle[6,],HAT0=masterData.Fsle[7,],GS0=masterData.Fsle[8,],
#                      BP0=masterData.Fsle[9,],BS0=masterData.Fsle[10,],JAX0=masterData.Fsle[11,])
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
# q25 = quantile(stack(FSLE)$values,probs=0.25)
# q75 = quantile(stack(FSLE)$values,probs=0.75)
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
# library(R.matlab)
# HYCOM = readMat('J:/Chpt_3/HYCOM_VerticalProfiles_forR.mat')
# chopInd = which(floor_date(as_date((HYCOM$Time[[1]][[1]]-719529),origin="1970-01-01"),unit="day")>as_date('2019-04-30',format="%Y-%m-%d"))[[1]]

# SSH ---------------------------------
# SSH = data.frame(HZ0=HYCOM$SSH[[1]][[1]][1:chopInd],OC0=HYCOM$SSH[[2]][[1]][1:chopInd],NC0=HYCOM$SSH[[3]][[1]][1:chopInd],
#                  BC0=HYCOM$SSH[[4]][[1]][1:chopInd],WC0=HYCOM$SSH[[5]][[1]][1:chopInd],NFC0=HYCOM$SSH[[6]][[1]][1:chopInd],
#                  HAT0=HYCOM$SSH[[7]][[1]][1:chopInd],GS0=HYCOM$SSH[[8]][[1]][1:chopInd],BP0=HYCOM$SSH[[9]][[1]][1:chopInd],
#                  BS0=HYCOM$SSH[[10]][[1]][1:chopInd],JAX0=HYCOM$SSH[[11]][[1]][1:chopInd])
# sites = colnames(SSH)
# 
# 
# # Histograms
# for (i in 1:11){
#   eval(parse(text=(paste(sites[i],' = ggplot(data=SSH)+geom_histogram(aes(x=',sites[i],'))+labs(x="m",y="",title=sites[i])',sep=""))))
# }
# png(file=paste(outDir,'/',"SSH_hist.png",sep=""),width = 800, height = 800, units = "px")
# grid.arrange(HZ0,OC0,NC0,BC0,WC0,NFC0,HAT0,GS0,BP0,BS0,JAX0, ncol=4,nrow=3,top="SSH")
# while (dev.cur()>1) {dev.off()}
# 
# # Plot boxplots and calculate quantiles to identify outliers
# ggplot(stack(SSH))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y="m",title="SSH")
# ggsave(paste(outDir,'/',"SSH_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)
# 
# q25 = quantile(stack(SSH)$values,probs=0.25)
# q75 = quantile(stack(SSH)$values,probs=0.75)
# iqr = q75-q25
# 
# # Remove outliers?
# for (i in 1:11){
#   q = round(SSH[,i],digits=10)
#   SSH[q==0,i] = NA
# }
# 
# # Check for missing dates
# missDate = list()
# for (i in 1:11){
# timeDiff = diff(HYCOM$Time[[i]][[1]][,])
# skipped = which(timeDiff>1)
# if(!is_empty(skipped)){
# missDate = c(missDate,list(skipped))
# } else {missDate = c(missDate,"No skipped dates")}
# }
# 
# # Interpolate missing dates
# fullSSH = data.frame(Time=seq.Date(from=as.Date("2016-02-01",origin="1970-01-01"),to=as.Date("2019-04-30",origin="1970-01-01"),by=1))
# for (i in 1:11){
#   fixedHYCOMTime = floor_date(as_date((HYCOM$Time[[i]][[1]][1:chopInd]-719529),origin="1970-01-01"),unit="day")
#   putWhere1 = match(fixedHYCOMTime,fullSSH$Time)
#   eval(parse(text=paste('fullSSH$',sites[i],'=NA',sep="")))
#   eval(parse(text=paste('fullSSH$',sites[i],'[putWhere1[!is.na(putWhere1)]]=SSH$',sites[i],'[-which(is.na(putWhere1))]',sep="")))
#   missSSH = apply(cbind(SSH[unlist(missDate[i]),i],SSH[unlist(missDate[i])+1,i]),MARGIN=1,mean)
#   missTimes = fixedHYCOMTime[unlist(missDate[i])]+1
#   putWhere2 = match(missTimes,fullSSH$Time) # indices where data is missing
#   missSSH = missSSH[-which(is.na(putWhere2))]
#   putWhere2 = putWhere2[-which(is.na(putWhere2))]
#   eval(parse(text=paste('fullSSH$',sites[i],'[putWhere2] = missSSH',sep="")))
# }
# 
# # Check for NAs and interpolate
# for (i in 1:11){
#   skippedBins = which(is.na(fullSSH[,i])) # missing data
#   missSSH = apply(cbind(fullSSH[skippedBins-1,i],fullSSH[skippedBins+1,i]),MARGIN=1,mean)
#   fullSSH[skippedBins,i] = missSSH
# }
# 
# # Create time lagged vectors
# startInd = which(fullSSH$Time==as.Date('2016-05-01',origin="1970-01-01"))
# fullLength = length(fullSSH$Time)
# for (i in 1:11){
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
#     eval(parse(text=paste('fullSSH$',sites[i],'Lag',lag,' = NA',sep="")))
#     eval(parse(text=paste('fullSSH$',sites[i],'Lag',lag,'[startInd:',fullLength,'] = fullSSH$',sites[i],'[lagInd:(',fullLength,'-lag)]',sep="")))
#   }
# }
# write.csv(fullSSH,file=paste(outDir,'/','SSH_TS.csv',sep=""),row.names=FALSE)



# Multi-depth Variables----------------
# vars = c('Salinity','Temperature','WaterU','WaterV')
# outDFs = c('SalFull','TempFull','UFull','VFull')
# depths = c(0,50,200,500,600,700,1000) # desired depth layers
# depInd = c(1,15,23,28,30,31,33)
# units = c("PPT",paste(intToUtf8(176),"C",sep=""),"m/s",'m/s')
# 
# for (j in 1:length(vars)){
#   
#   eval(parse(text=paste('full',vars[j],' = data.frame(Time=seq.Date(from=as.Date("2016-02-01",origin="1970-01-01"),to=as.Date("2019-04-30",origin="1970-01-01"),by=1))',sep="")))
#   
#   for (l in 1:length(depths)){
#     
#     eval(parse(text=paste('Temp = data.frame(HZ',depths[l],'=HYCOM$',vars[j],'[[1]][[1]][',depInd[l],
#                           ',1:chopInd],OC',depths[l],'=HYCOM$',vars[j],'[[2]][[1]][',depInd[l],
#                           ',1:chopInd],NC',depths[l],'=HYCOM$',vars[j],'[[3]][[1]][',depInd[l],
#                           ',1:chopInd],BC',depths[l],'=HYCOM$',vars[j],'[[4]][[1]][',depInd[l],
#                           ',1:chopInd],WC',depths[l],'=HYCOM$',vars[j],'[[5]][[1]][',depInd[l],
#                           ',1:chopInd],NFC',depths[l],'=HYCOM$',vars[j],'[[6]][[1]][',depInd[l],
#                           ',1:chopInd],HAT',depths[l],'=HYCOM$',vars[j],'[[7]][[1]][',depInd[l],
#                           ',1:chopInd],GS',depths[l],'=HYCOM$',vars[j],'[[8]][[1]][',depInd[l],
#                           ',1:chopInd],BP',depths[l],'=HYCOM$',vars[j],'[[9]][[1]][',depInd[l],
#                           ',1:chopInd],BS',depths[l],'=HYCOM$',vars[j],'[[10]][[1]][',depInd[l],
#                           ',1:chopInd],JAX',depths[l],'=HYCOM$',vars[j],'[[11]][[1]][',depInd[l],
#                           ',1:chopInd])',sep="")))
#     sites = colnames(Temp)
#     siteName = str_replace(sites,as.character(depths[l]),"")
#     
#     for (i in 1:11){
#       eval(parse(text=paste('Temp$',sites[i],'[is.nan(Temp$',sites[i], ')] = NA',sep="")))
#       if (vars[j]=="Temperature"){
#         eval(parse(text=paste('Temp$',sites[i],'[Temp$',sites[i],'==0] = NA',sep="")))
#       }
#       if (vars[j]=="WaterU"){
#         eval(parse(text=paste('Temp$',sites[i],'[Temp$',sites[i],'==0] = NA',sep="")))
#       }
#       if (vars[j]=="WaterU"){
#         eval(parse(text=paste('Temp$',sites[i],'[Temp$',sites[i],'==0] = NA',sep="")))
#       }
#     }
#     
#     # Histograms
#     for (i in 1:11){
#       eval(parse(text=(paste(siteName[i],' = ggplot(data=Temp)+geom_histogram(aes(x=',
#                              sites[i],'))+labs(x=units[j],title=sites[i])',sep=""))))
#     }
#     png(file=paste(outDir,'/',vars[j],'_',as.character(depths[l]),"_hist.png",sep=""),width = 800, height = 800, units = "px")
#     grid.arrange(HZ,NC,OC,BC,WC,NFC,HAT,GS,BP,BS,JAX, ncol=4,nrow=3,top=paste(vars[j],' at ',as.character(depths[l]),'m',sep=""))
#     while (dev.cur()>1) {dev.off()}
#     
#     # Plot boxplots and calculate quantiles to identify outliers
#     ggplot(stack(Temp))+geom_boxplot(aes(x=ind,y=values))+labs(x="Site & Depth",y="m",title=paste(vars[j]),' at ',as.character(depths[l]),'m',sep="")
#     ggsave(paste(outDir,'/',vars[j],'_',as.character(depths[l]),"_boxplot.png",sep=""),device="png",width=600,height=800,units="px",scale=5,dpi=600)
#     
#     # Remove outliers?
#     q25 = quantile(stack(Temp[,1:10])$values,probs=0.25,na.rm=TRUE)
#     q75 = quantile(stack(Temp[,1:10])$values,probs=0.75,na.rm=TRUE)
#     iqr = q75-q25
#     upperThresh = q75+(1.5*iqr)
#     lowerThresh = q25-(1.5*iqr)
#     
#     # Check for missing dates
#     missDate = list()
#     for (i in 1:11){
#       timeDiff = diff(HYCOM$Time[[i]][[1]][,])
#       skipped = which(timeDiff>1)
#       if(!is_empty(skipped)){
#         missDate = c(missDate,list(skipped))
#       } else {missDate = c(missDate,"No skipped dates")}
#     }
#     
#     # Interpolate missing dates
#     
#     for (i in 1:11){
#       fixedHYCOMTime = floor_date(as_date((HYCOM$Time[[i]][[1]][1:chopInd]-719529),origin="1970-01-01"),unit="day")
#       eval(parse(text=paste('putWhere1 = match(fixedHYCOMTime,full',vars[j],'$Time)',sep="")))
#       eval(parse(text=paste('full',vars[j],'$',sites[i],'=NA',sep="")))
#       eval(parse(text=paste('full',vars[j],'$',sites[i],'[putWhere1[!is.na(putWhere1)]]=Temp$',sites[i],'[-which(is.na(putWhere1))]',sep="")))
#       missData = apply(cbind(Temp[unlist(missDate[i]),i],Temp[unlist(missDate[i])+1,i]),MARGIN=1,mean)
#       missTimes = fixedHYCOMTime[unlist(missDate[i])]+1
#       eval(parse(text=paste('putWhere2 = match(missTimes,full',vars[j],'$Time)',sep=""))) # indices where data is missing
#       missData = missData[-which(is.na(putWhere2))]
#       putWhere2 = putWhere2[-which(is.na(putWhere2))]
#       eval(parse(text=paste('full',vars[j],'$',sites[i],'[putWhere2] = missData',sep="")))
#       
#     }
#   }
#   
#   # Check for NAs and interpolate (Need a catch for 2+ consecutive data points missing)
#   for (i in 2:78){
#     eval(parse(text=paste('skippedBins = which(is.na(full',vars[j],'[,i]))',sep=""))) # missing data
#     if(any(diff(skippedBins)==1)){
#       return("WARNING: 2 OR MORE NAs IN A ROW; CAN'T INTERPOLATE")
#     }
#     eval(parse(text=paste('missData = apply(cbind(full',vars[j],'[skippedBins-1,i],full',vars[j],'[skippedBins+1,i]),MARGIN=1,mean)',sep="")))
#     eval(parse(text=paste('full',vars[j],'[skippedBins,i] = missData',sep="")))
#   }
# 
#   # Create time lagged vectors
#   eval(parse(text=paste('startInd = which(full',vars[j],'$Time==as.Date("2016-05-01",origin="1970-01-01"))',sep="")))
#   eval(parse(text=paste('fullLength = length(full',vars[j],'$Time)',sep="")))
#   eval(parse(text=paste('fullColNames = colnames(full',vars[j],')',sep="")))
#   for (i in 2:length(fullColNames)){
#     for (k in 1:3){
#       if (k==1){
#         lag = 7
#         lagInd = startInd-lag
#       } else if (k==2){
#         lag=14
#         lagInd = startInd-lag
#       } else if (k==3){
#         lag=21
#         lagInd = startInd-lag
#       }
#       # baseData = paste(sites[i],'0',sep="")
#       # thisDepth = 
#       eval(parse(text=paste('full',vars[j],'$',fullColNames[i],'Lag',lag,' = NA',sep="")))
#       eval(parse(text=paste('full',vars[j],'$',fullColNames[i],'Lag',lag,'[startInd:',fullLength,'] = full',vars[j],'$',fullColNames[i],'[lagInd:(',fullLength,'-lag)]',sep="")))
#     }
#   }
#   
#   # Save data as .csv
#   eval(parse(text=paste('write.csv(full',vars[j],',file=',outDir,'/',vars[j],'_TS.csv)',sep="")),row.names=FALSE)
#   
# }

## Clean Chla Data  --------------------------------------
load('J:/Chpt_3/Chl_TS.Rdata')
Chl = data.frame(Time=as.Date(masterData.Time[1,],origin="1970-01-01"),HZ0=masterData.Chl[1,],OC0=masterData.Chl[2,],NC0=masterData.Chl[3,],BC0=masterData.Chl[4,],
                  WC0=masterData.Chl[5,],NFC0=masterData.Chl[6,],HAT0=masterData.Chl[7,],GS0=masterData.Chl[8,],
                  BP0=masterData.Chl[9,],BS0=masterData.Chl[10,],JAX0=masterData.Chl[11,])

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
  missChlLin = approx(x=Chl[datBins,1],y=Chl[datBins,i],xout=Chl[,1],rule=1,method="linear",ties=mean)

  ChlBasis = bSpline(as.numeric(Chl[,1]),df=20,Boundary.knots=c(min(Chl[,1],na.rm=TRUE),max(Chl[,1],na.rm=TRUE)))
  ChlMod = gam(Chl[,i]~s(as.numeric(Chl[,1]),bs="cs",k=21),family=gaussian)
  missChlBS = ChlBasis%*%coef(ChlMod)[2:21]
  
  missChlSS = smooth.spline(x=Chl[datBins,1],y=Chl[datBins,i],all.knots=TRUE)
  
  missChlNS = spline(x=Chl[datBins,1],y=Chl[datBins,i],xout=Chl[,1],method="natural")
  
  Chl[,i] = missChl
}

# Create time lagged vectors

## Make modeling csv files ---------------------------

