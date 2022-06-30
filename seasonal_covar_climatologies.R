library(stringr)
library(lubridate)
library(EFDR)

outDir = 'J:/Chpt_3/Predictions'
predictWinter = list()
predictSpring = list()
predictSummer = list()
predictFall = list()

# HYCOM vars -----------------------------------
varList = c("EKE","SSH","Salinity","Temperature","VelocityAsp","VelocityMag")
depths = c('_0_','_200_','_400_','_700_')
inDir = 'J:/Chpt_3/HYCOM/0.08deg'
fileList = dir(inDir,".Rdata",full.names=TRUE,recursive=FALSE)
latMat = numeric()
lonMat = numeric()

for (i in 1:length(varList)){
  
  if (varList[i]=='SSH' | varList[i]=="EKE"){
    depth = '_0_'
  } else {
    depth = depths
  }
  
  for (j in 1:length(depth)){
    
    varFiles = which(str_detect(fileList,varList[i]) & str_detect(fileList,depth[j]))
    eval(parse(text=paste('thisVar_Winter = numeric()',sep="")))
    eval(parse(text=paste('thisVar_Spring = numeric()',sep="")))
    eval(parse(text=paste('thisVar_Summer = numeric()',sep="")))
    eval(parse(text=paste('thisVar_Fall = numeric()',sep="")))
    varName = str_remove_all(paste(varList[i],depth[j],sep=""),'_')
    
    
    for (k in 1:length(varFiles)){
      # get 6-digit datestamps from file names
      fileDate = str_extract(fileList[varFiles[k]],"\\d\\d\\d\\d\\d\\d\\d\\d") 
      time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                        str_sub(fileDate,start=5L,end=6L),'-',
                        str_sub(fileDate,start=7L,end=8L),sep="")
      
      thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
      
      if (month(thisTime)==1){
        load(fileList[varFiles[k]])
        dataVec = stack(data.frame(data))[,1]
        eval(parse(text=paste('thisVar_Winter = rbind(thisVar_Winter,dataVec)',sep="")))
        latMat = cbind(latMat,lats)
        lonMat = cbind(lonMat,lons)
      } else if (month(thisTime)==4){
        load(fileList[varFiles[k]])
        dataVec = stack(data.frame(data))[,1]
        eval(parse(text=paste('thisVar_Spring = rbind(thisVar_Spring,dataVec)',sep="")))
        latMat = cbind(latMat,lats)
        lonMat = cbind(lonMat,lons)
      } else if (month(thisTime)==7){
        load(fileList[varFiles[k]])
        dataVec = stack(data.frame(data))[,1]
        eval(parse(text=paste('thisVar_Summer = rbind(thisVar_Summer,dataVec)',sep="")))
        latMat = cbind(latMat,lats)
        lonMat = cbind(lonMat,lons)
      } else if (month(thisTime)==10){
        load(fileList[varFiles[k]])
        dataVec = stack(data.frame(data))[,1]
        eval(parse(text=paste('thisVar_Fall = rbind(thisVar_Fall,dataVec)',sep="")))
        latMat = cbind(latMat,lats)
        lonMat = cbind(lonMat,lons)
      }
      
    }
    
    matplot(latMat,main=paste(varName," Latitudes"),type="p")
    matplot(lonMat,main=paste(varName," Longitudes"),type="p")
    
    predictWinter$lat = rep(lats,length.out=length(lons)*length(lats))
    predictWinter$lon = rep(lons,each=length(lats))
    eval(parse(text=paste('predictWinter$',varName,'=apply(thisVar_Winter,MARGIN=2,mean,na.rm=TRUE)',sep="")))
    predictSpring$lat = rep(lats,length.out=length(lons)*length(lats))
    predictSpring$lon = rep(lons,each=length(lats))
    eval(parse(text=paste('predictSpring$',varName,'=apply(thisVar_Spring,MARGIN=2,mean,na.rm=TRUE)',sep="")))
    predictSummer$lat = rep(lats,length.out=length(lons)*length(lats))
    predictSummer$lon = rep(lons,each=length(lats))
    eval(parse(text=paste('predictSummer$',varName,'=apply(thisVar_Summer,MARGIN=2,mean,na.rm=TRUE)',sep="")))
    predictFall$lat = rep(lats,length.out=length(lons)*length(lats))
    predictFall$lon = rep(lons,each=length(lats))
    eval(parse(text=paste('predictFall$',varName,'=apply(thisVar_Fall,MARGIN=2,mean,na.rm=TRUE)',sep="")))
    
  }
}

save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))

# FSLE ------------------------------------

inDir = 'J:/Chpt_3/FSLE/0.08deg'
fileList = dir(inDir,".Rdata",full.names=TRUE,recursive=FALSE)

FSLE_Winter = numeric()
FSLE_Spring = numeric()
FSLE_Summer = numeric()
FSLE_Fall = numeric()

for (k in 1:length(fileList)-1){
  # get 6-digit datestamps from file names
  fileDate = str_extract(fileList[k],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  if (month(thisTime)==1){
    load(fileList[k])
    dataVec = stack(data.frame(data))[,1]
    FSLE_Winter = rbind(FSLE_Winter,dataVec)
  } else if (month(thisTime)==4){
    load(fileList[k])
    dataVec = stack(data.frame(data))[,1]
    FSLE_Spring = rbind(FSLE_Spring,dataVec)  
    } else if (month(thisTime)==7){
    load(fileList[k])
    dataVec = stack(data.frame(data))[,1]
    FSLE_Summer = rbind(FSLE_Summer,dataVec)
  } else if (month(thisTime)==10){
    load(fileList[k])
    dataVec = stack(data.frame(data))[,1]
    FSLE_Fall = rbind(FSLE_Fall,dataVec)
  }
  
}

predictWinter$FSLE0 = apply(FSLE_Winter,2,mean,na.rm=TRUE)
predictSpring$FSLE0 = apply(FSLE_Spring,2,mean,na.rm=TRUE)
predictSummer$FSLE0 = apply(FSLE_Summer,2,mean,na.rm=TRUE)
predictFall$FSLE0 = apply(FSLE_Fall,2,mean,na.rm=TRUE)

save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))

# Chl ------------------------------------------------------

inDir = 'J:/Chpt_3/Chla/0.0466deg'
fileList = dir(inDir,".Rdata",full.names=TRUE,recursive=FALSE)

Chl_Winter = numeric()
Chl_Spring = numeric()
Chl_Summer = numeric()
Chl_Fall = numeric()

for (k in 4:length(fileList)){
  # get 6-digit datestamps from file names
  fileDate = str_extract(fileList[k],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  if (month(thisTime)==1){
    load(fileList[k])
    dataVec = stack(data.frame(data))[,1]
    Chl_Winter = rbind(Chl_Winter,dataVec)
  } else if (month(thisTime)==4){
    load(fileList[k])
    dataVec = stack(data.frame(data))[,1]
    Chl_Spring = rbind(Chl_Spring,dataVec)  
  } else if (month(thisTime)==7){
    load(fileList[k])
    dataVec = stack(data.frame(data))[,1]
    Chl_Summer = rbind(Chl_Summer,dataVec)
  } else if (month(thisTime)==10){
    load(fileList[k])
    dataVec = stack(data.frame(data))[,1]
    Chl_Fall = rbind(Chl_Fall,dataVec)
  }
  
}


# regrid to 0.08deg grid
newx = seq(278,297,by=0.08)
newy = seq(24,44,by=0.08)
ChlWinterDF = data.frame(z=stack(data.frame(apply(Chl_Winter,2,mean,na.rm=TRUE)))[,1],
                         y=rep(lats,length.out=length(lons)*length(lats)),
                         x=rep(lons,each=length(lats)))
ChlSpringDF = data.frame(z=stack(data.frame(apply(Chl_Spring,2,mean,na.rm=TRUE)))[,1],
                         y=rep(lats,length.out=length(lons)*length(lats)),
                         x=rep(lons,each=length(lats)))
ChlSummerDF = data.frame(z=stack(data.frame(apply(Chl_Summer,2,mean,na.rm=TRUE)))[,1],
                         y=rep(lats,length.out=length(lons)*length(lats)),
                         x=rep(lons,each=length(lats)))
ChlFallDF = data.frame(z=stack(data.frame(apply(Chl_Fall,2,mean,na.rm=TRUE)))[,1],
                         y=rep(lats,length.out=length(lons)*length(lats)),
                         x=rep(lons,each=length(lats)))

ChlWinterDF$z[is.na(ChlWinterDF$z)] = 0
newWin = regrid(ChlWinterDF,n1=length(newx),n2=length(newy),method="idw")
ChlSpringDF$z[is.na(ChlSpringDF$z)] = 0
newSpr = regrid(ChlSpringDF,n1=length(newx),n2=length(newy),method="idw")
ChlSummerDF$z[is.na(ChlSummerDF$z)] = 0
newSum = regrid(ChlSummerDF,n1=length(newx),n2=length(newy),method="idw")
ChlFallDF$z[is.na(ChlFallDF$z)] = 0
newFall = regrid(ChlFallDF,n1=length(newx),n2=length(newy),method="idw")

predictWinter$Chl0 = unlist(stack(data.frame(matrix(newWin$z,ncol=length(newx),byrow=TRUE)))[,1])
predictSpring$Chl0 = unlist(stack(data.frame(matrix(newSpr$z,ncol=length(newx),byrow=TRUE)))[,1])
predictSummer$Chl0 = unlist(stack(data.frame(matrix(newSum$z,ncol=length(newx),byrow=TRUE)))[,1])
predictFall$Chl0 = unlist(stack(data.frame(matrix(newFall$z,ncol=length(newx),byrow=TRUE)))[,1])

save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))

# Eddy Dist -----------------------------------------------------

inDir = 'J:/Chpt_3/Eddies/Grids'

# anticyclonic eddies
AfileList = dir(inDir,"AEddy",full.names=TRUE,recursive=FALSE)
AEddy_Winter = numeric()
AEddy_Spring = numeric()
AEddy_Summer = numeric()
AEddy_Fall = numeric()

for (k in 1:length(AfileList)){
  # get 6-digit datestamps from file names
  fileDate = str_extract(AfileList[k],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  if (month(thisTime)==1){
    load(AfileList[k])
    dataVec = stack(data.frame(data))[,1]
    AEddy_Winter = rbind(AEddy_Winter,dataVec)
  } else if (month(thisTime)==4){
    load(AfileList[k])
    dataVec = stack(data.frame(data))[,1]
    AEddy_Spring = rbind(AEddy_Spring,dataVec)  
  } else if (month(thisTime)==7){
    load(AfileList[k])
    dataVec = stack(data.frame(data))[,1]
    AEddy_Summer = rbind(AEddy_Summer,dataVec)
  } else if (month(thisTime)==10){
    load(AfileList[k])
    dataVec = stack(data.frame(data))[,1]
    AEddy_Fall = rbind(AEddy_Fall,dataVec)
  }
  
}

predictWinter$AEddyDist0 = apply(AEddy_Winter,2,mean,na.rm=TRUE)
predictSpring$AEddyDist0 = apply(AEddy_Spring,2,mean,na.rm=TRUE)
predictSummer$AEddyDist0 = apply(AEddy_Summer,2,mean,na.rm=TRUE)
predictFall$AEddyDist0 = apply(AEddy_Fall,2,mean,na.rm=TRUE)

# Cyclonic eddies
CfileList = dir(inDir,"CEddy",full.names=TRUE,recursive=FALSE)

CEddy_Winter = numeric()
CEddy_Spring = numeric()
CEddy_Summer = numeric()
CEddy_Fall = numeric()

for (k in 1:length(CfileList)){
  # get 6-digit datestamps from file names
  fileDate = str_extract(CfileList[k],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  time_temp = paste(str_sub(fileDate,start=1L,end=4L),'-',
                    str_sub(fileDate,start=5L,end=6L),'-',
                    str_sub(fileDate,start=7L,end=8L),sep="")
  
  thisTime = as.Date(time_temp,format='%Y-%m-%d',tz="UTC")
  
  if (month(thisTime)==1){
    load(CfileList[k])
    dataVec = stack(data.frame(data))[,1]
    CEddy_Winter = rbind(CEddy_Winter,dataVec)
  } else if (month(thisTime)==4){
    load(CfileList[k])
    dataVec = stack(data.frame(data))[,1]
    CEddy_Spring = rbind(CEddy_Spring,dataVec)  
  } else if (month(thisTime)==7){
    load(CfileList[k])
    dataVec = stack(data.frame(data))[,1]
    CEddy_Summer = rbind(CEddy_Summer,dataVec)
  } else if (month(thisTime)==10){
    load(CfileList[k])
    dataVec = stack(data.frame(data))[,1]
    CEddy_Fall = rbind(CEddy_Fall,dataVec)
  }
  
}

predictWinter$CEddyDist0 = apply(CEddy_Winter,2,mean,na.rm=TRUE)
predictSpring$CEddyDist0 = apply(CEddy_Spring,2,mean,na.rm=TRUE)
predictSummer$CEddyDist0 = apply(CEddy_Summer,2,mean,na.rm=TRUE)
predictFall$CEddyDist0 = apply(CEddy_Fall,2,mean,na.rm=TRUE)

save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))

# Convert prediction data to data frames -----------------

predictWinter = data.frame(predictWinter)
predictSpring = data.frame(predictSpring)
predictSummer = data.frame(predictSummer)
predictFall = data.frame(predictFall)

# rename columns to match model covars
colnames(predictWinter) = c("lat","lon","EKE0","SSH0","Sal0","Sal200","Sal400","Sal700",
                            "Temp0","Temp200","Temp400","Temp700","VelAsp0","VelAsp200",
                            "VelAsp400","VelAsp700","VelMag0","VelMag200","VelMag400",
                            "VelMag700","FSLE0","Chl0","AEddyDist0","CEddyDist0")
colnames(predictSpring) = c("lat","lon","EKE0","SSH0","Sal0","Sal200","Sal400","Sal700",
                            "Temp0","Temp200","Temp400","Temp700","VelAsp0","VelAsp200",
                            "VelAsp400","VelAsp700","VelMag0","VelMag200","VelMag400",
                            "VelMag700","FSLE0","Chl0","AEddyDist0","CEddyDist0")
colnames(predictSummer) = c("lat","lon","EKE0","SSH0","Sal0","Sal200","Sal400","Sal700",
                            "Temp0","Temp200","Temp400","Temp700","VelAsp0","VelAsp200",
                            "VelAsp400","VelAsp700","VelMag0","VelMag200","VelMag400",
                            "VelMag700","FSLE0","Chl0","AEddyDist0","CEddyDist0")
colnames(predictFall) = c("lat","lon","EKE0","SSH0","Sal0","Sal200","Sal400","Sal700",
                            "Temp0","Temp200","Temp400","Temp700","VelAsp0","VelAsp200",
                            "VelAsp400","VelAsp700","VelMag0","VelMag200","VelMag400",
                            "VelMag700","FSLE0","Chl0","AEddyDist0","CEddyDist0")

# transform data as necessary
# -Inf will occur in Chl0 + FSLE0 due to 0 values, set to NA

predictWinter$log_Chl0 = log10(predictWinter$Chl0)
predictWinter$log_Chl0[is.infinite(predictWinter$log_Chl0)] = NA
predictWinter$log_abs_FSLE0 = log10(abs(predictWinter$FSLE0))
predictWinter$log_abs_FSLE0[is.infinite(predictWinter$log_abs_FSLE0)] = NA
predictWinter$sqrt_CEddyDist0 = sqrt(predictWinter$CEddyDist0)
predictWinter$sqrt_AEddyDist0 = sqrt(predictWinter$AEddyDist0)
predictWinter$sqrt_VelAsp0 = sqrt(predictWinter$VelAsp0)
predictWinter$sqrt_VelAsp700 = sqrt(predictWinter$VelAsp700)

predictSpring$log_Chl0 = log10(predictSpring$Chl0)
predictSpring$log_Chl0[is.infinite(predictSpring$log_Chl0)] = NA
predictSpring$log_abs_FSLE0 = log10(abs(predictSpring$FSLE0))
predictSpring$log_abs_FSLE0[is.infinite(predictSpring$log_abs_FSLE0)] = NA
predictSpring$sqrt_CEddyDist0 = sqrt(predictSpring$CEddyDist0)
predictSpring$sqrt_AEddyDist0 = sqrt(predictSpring$AEddyDist0)
predictSpring$sqrt_VelAsp0 = sqrt(predictSpring$VelAsp0)
predictSpring$sqrt_VelAsp700 = sqrt(predictSpring$VelAsp700)

predictSummer$log_Chl0 = log10(predictSummer$Chl0)
predictSummer$log_Chl0[is.infinite(predictSummer$log_Chl0)] = NA
predictSummer$log_abs_FSLE0 = log10(abs(predictSummer$FSLE0))
predictSummer$log_abs_FSLE0[is.infinite(predictSummer$log_abs_FSLE0)] = NA
predictSummer$sqrt_CEddyDist0 = sqrt(predictSummer$CEddyDist0)
predictSummer$sqrt_AEddyDist0 = sqrt(predictSummer$AEddyDist0)
predictSummer$sqrt_VelAsp0 = sqrt(predictSummer$VelAsp0)
predictSummer$sqrt_VelAsp700 = sqrt(predictSummer$VelAsp700)

predictFall$log_Chl0 = log10(predictFall$Chl0)
predictFall$log_Chl0[is.infinite(predictFall$log_Chl0)] = NA
predictFall$log_abs_FSLE0 = log10(abs(predictFall$FSLE0))
predictFall$log_abs_FSLE0[is.infinite(predictFall$log_abs_FSLE0)] = NA
predictFall$sqrt_CEddyDist0 = sqrt(predictFall$CEddyDist0)
predictFall$sqrt_AEddyDist0 = sqrt(predictFall$AEddyDist0)
predictFall$sqrt_VelAsp0 = sqrt(predictFall$VelAsp0)
predictFall$sqrt_VelAsp700 = sqrt(predictFall$VelAsp700)

# get elevation values to later remove prediction values on land
load('J:/Chpt_3/GEBCO/DownsampledGrid_08deg.Rdata')
x = seq(278,297,by=0.08)
y = seq(24,44,by=0.08)
depth = unlist(stack(data.frame(matrix(newGrid$Depth,ncol=length(x),byrow=TRUE)))[,1])
depth[depth>=0] = NA
predictWinter$Depth = depth
predictSpring$Depth = depth
predictSummer$Depth = depth
predictFall$Depth = depth

save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))
