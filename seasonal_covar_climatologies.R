library(stringr)
library(lubridate)
library(EFDR)

outDir = 'E:/Chpt_3/Predictions'
# predictWinter = list()
# predictSpring = list()
# predictSummer = list()
# predictFall = list()

# HYCOM vars -----------------------------------
varList = c("EKE","SSH","Salinity","Temperature","VelocityAsp","VelocityMag")
# depths = c('_0_','_200_','_400_','_700_')
depths = '_0_'
inDir = 'E:/Chpt_3/HYCOM/0.08deg_Surface'
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
    # thisVar_Winter = numeric()
    # thisVar_Spring = numeric()
    # thisVar_Summer = numeric()
    # thisVar_Fall = numeric()
    
    varName = str_remove_all(paste(varList[i],depth[j],sep=""),'_')
    
    for (k in 1:length(varFiles)){
      # get 6-digit datestamps from file names
      fileDate = str_extract(fileList[varFiles[k]],"\\d{8}")
      thisTime = as_datetime(fileDate,format='%Y%m%d',tz="UTC")
      thisMonth = month(thisTime)
      
      # if (month(thisTime)==1){
      #   load(fileList[varFiles[k]])
      #   dataVec = stack(data.frame(data))[,1]
      #   thisVar_Winter = rbind(thisVar_Winter,dataVec)
      #   latMat = cbind(latMat,lats)
      #   lonMat = cbind(lonMat,lons)
      # } else if (month(thisTime)==4){
      #   load(fileList[varFiles[k]])
      #   dataVec = stack(data.frame(data))[,1]
      #   thisVar_Spring = rbind(thisVar_Spring,dataVec)
      #   latMat = cbind(latMat,lats)
      #   lonMat = cbind(lonMat,lons)
      # } else if (month(thisTime)==7){
      #   load(fileList[varFiles[k]])
      #   dataVec = stack(data.frame(data))[,1]
      #   thisVar_Summer = rbind(thisVar_Summer,dataVec)
      #   latMat = cbind(latMat,lats)
      #   lonMat = cbind(lonMat,lons)
      # } else if (month(thisTime)==10){
      #   load(fileList[varFiles[k]])
      #   dataVec = stack(data.frame(data))[,1]
      #   thisVar_Fall = rbind(thisVar_Fall,dataVec)
      #   latMat = cbind(latMat,lats)
      #   lonMat = cbind(lonMat,lons)
      # }
      
      load(fileList[varFiles[k]])
      dataVec = stack(data.frame(data))[,1]
      if (!exists(eval(parse(text=paste('\'thisVar_',thisMonth,'\'',sep=""))))){
        eval(parse(text=paste('thisVar_',thisMonth,' = dataVec',sep="")))
      }else{
        eval(parse(text=paste('thisVar_',thisMonth,' = rbind(thisVar_',thisMonth,',dataVec)',sep="")))
      }
      
    }
    
    # predictWinter$lat = rep(lats,length.out=length(lons)*length(lats))
    # predictWinter$lon = rep(lons,each=length(lats))
    # eval(parse(text=paste('predictWinter$',varName,'=apply(thisVar_Winter,MARGIN=2,mean,na.rm=TRUE)',sep="")))
    # predictSpring$lat = rep(lats,length.out=length(lons)*length(lats))
    # predictSpring$lon = rep(lons,each=length(lats))
    # eval(parse(text=paste('predictSpring$',varName,'=apply(thisVar_Spring,MARGIN=2,mean,na.rm=TRUE)',sep="")))
    # predictSummer$lat = rep(lats,length.out=length(lons)*length(lats))
    # predictSummer$lon = rep(lons,each=length(lats))
    # eval(parse(text=paste('predictSummer$',varName,'=apply(thisVar_Summer,MARGIN=2,mean,na.rm=TRUE)',sep="")))
    # predictFall$lat = rep(lats,length.out=length(lons)*length(lats))
    # predictFall$lon = rep(lons,each=length(lats))
    # eval(parse(text=paste('predictFall$',varName,'=apply(thisVar_Fall,MARGIN=2,mean,na.rm=TRUE)',sep="")))
    
    for (l in 1:12){
      if(i==1){
        eval(parse(text=paste('predict_',l,' = data.frame(lat = rep(lats,length.out=length(lons)*length(lats)))',sep="")))
        eval(parse(text=paste('predict_',l,'$lon = rep(lons,each=length(lats))',sep="")))
      }
      eval(parse(text=paste('predict_',l,'$',varName,' = apply(thisVar_',l,',MARGIN=2,mean,na.rm=TRUE)',sep="")))
    }
    
  }
}

# save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/SeasonalPredictionData.Rdata',sep=""))
save(predict_1,predict_2,predict_3,predict_4,predict_5,predict_6,
     predict_7,predict_8,predict_9,predict_10,predict_11,predict_12,
     file=paste(outDir,'/MonthlyPredictionData.Rdata',sep=""))

# FSLE ------------------------------------

inDir = 'E:/Chpt_3/FSLE/0.08deg'
fileList = dir(inDir,".Rdata",full.names=TRUE,recursive=FALSE)

# FSLE_Winter = numeric()
# FSLE_Spring = numeric()
# FSLE_Summer = numeric()
# FSLE_Fall = numeric()

for (k in 1:(length(fileList)-1)){
  # get 6-digit datestamps from file names
  fileDate = str_extract(fileList[k],'\\d{8}') 
  thisTime = as_datetime(fileDate,format='%Y%m%d',tz="UTC")
  thisMonth = month(thisTime)
  
  load(fileList[k])
  
  # if (month(thisTime)==1){
  #   load(fileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   FSLE_Winter = rbind(FSLE_Winter,dataVec)
  # } else if (month(thisTime)==4){
  #   load(fileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   FSLE_Spring = rbind(FSLE_Spring,dataVec)  
  #   } else if (month(thisTime)==7){
  #   load(fileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   FSLE_Summer = rbind(FSLE_Summer,dataVec)
  # } else if (month(thisTime)==10){
  #   load(fileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   FSLE_Fall = rbind(FSLE_Fall,dataVec)
  # }
  
  dataVec = stack(data.frame(data))[,1]
  if (!exists(eval(parse(text=paste('\'FSLE_',thisMonth,'\'',sep=""))))){
    eval(parse(text=paste('FSLE_',thisMonth,' = dataVec',sep="")))
  }else{
    eval(parse(text=paste('FSLE_',thisMonth,' = rbind(FSLE_',thisMonth,',dataVec)',sep="")))
  }
  
}

# predictWinter$FSLE0 = apply(FSLE_Winter,2,mean,na.rm=TRUE)
# predictSpring$FSLE0 = apply(FSLE_Spring,2,mean,na.rm=TRUE)
# predictSummer$FSLE0 = apply(FSLE_Summer,2,mean,na.rm=TRUE)
# predictFall$FSLE0 = apply(FSLE_Fall,2,mean,na.rm=TRUE)

for (l in 1:12){
  eval(parse(text=paste('predict_',l,' = data.frame(FSLE = apply(FSLE_',l,',MARGIN=2,mean,na.rm=TRUE))',sep="")))
}

# save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/SeasonalPredictionData.Rdata',sep=""))
save(predict_1,predict_2,predict_3,predict_4,predict_5,predict_6,
     predict_7,predict_8,predict_9,predict_10,predict_11,predict_12,
     file=paste(outDir,'/MonthlyPredictionData.Rdata',sep=""))

# Chl ------------------------------------------------------

inDir = 'E:/Chpt_3/Chla/0.0466deg'
fileList = dir(inDir,".Rdata",full.names=TRUE,recursive=FALSE)

# Chl_Winter = numeric()
# Chl_Spring = numeric()
# Chl_Summer = numeric()
# Chl_Fall = numeric()

for (k in 4:length(fileList)){
  # get 6-digit datestamps from file names
  fileDate = str_extract(fileList[k],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  thisTime = as_datetime(fileDate,format='%Y%m%d',tz="UTC")
  thisMonth = month(thisTime)
  
  # if (month(thisTime)==1){
  #   load(fileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   Chl_Winter = rbind(Chl_Winter,dataVec)
  # } else if (month(thisTime)==4){
  #   load(fileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   Chl_Spring = rbind(Chl_Spring,dataVec)  
  # } else if (month(thisTime)==7){
  #   load(fileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   Chl_Summer = rbind(Chl_Summer,dataVec)
  # } else if (month(thisTime)==10){
  #   load(fileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   Chl_Fall = rbind(Chl_Fall,dataVec)
  # }
  
  dataVec = stack(data.frame(data))[,1]
  if (!exists(eval(parse(text=paste('\'Chl_',thisMonth,'\'',sep=""))))){
    eval(parse(text=paste('Chl_',thisMonth,' = dataVec',sep="")))
  }else{
    eval(parse(text=paste('Chl_',thisMonth,' = rbind(Chl_',thisMonth,',dataVec)',sep="")))
  }
  
}


# regrid to 0.08deg grid
newx = seq(278,297,by=0.08)
newy = seq(24,44,by=0.08)
for (i in 1:12){
  eval(parse(text=paste('Chl',i,'DF = data.frame(z=stack(data.frame(apply(Chl_',i,',2,mean,na.rm=TRUE)))[,1],
                           y=rep(lats,length.out=length(lons)*length(lats)),
                           x=rep(lons,each=length(lats)))',sep="")))
  eval(parse(text=paste('ChlDF',i,'$z[is.na(ChlDF',i,'$z)] = 0',sep="")))
  eval(parse(paste('new_',i,' = regrid(ChlDF',i,',n1=length(newx),n2=length(newy),method="idw")',sep="")))
  paste('predict_',i,'$Chl0 = unlist(stack(data.frame(matrix(new_',i,'$z,ncol=length(newx),byrow=TRUE)))[,1])',sep="")
}

# ChlWinterDF = data.frame(z=stack(data.frame(apply(Chl_Winter,2,mean,na.rm=TRUE)))[,1],
#                          y=rep(lats,length.out=length(lons)*length(lats)),
#                          x=rep(lons,each=length(lats)))
# ChlSpringDF = data.frame(z=stack(data.frame(apply(Chl_Spring,2,mean,na.rm=TRUE)))[,1],
#                          y=rep(lats,length.out=length(lons)*length(lats)),
#                          x=rep(lons,each=length(lats)))
# ChlSummerDF = data.frame(z=stack(data.frame(apply(Chl_Summer,2,mean,na.rm=TRUE)))[,1],
#                          y=rep(lats,length.out=length(lons)*length(lats)),
#                          x=rep(lons,each=length(lats)))
# ChlFallDF = data.frame(z=stack(data.frame(apply(Chl_Fall,2,mean,na.rm=TRUE)))[,1],
#                          y=rep(lats,length.out=length(lons)*length(lats)),
#                          x=rep(lons,each=length(lats)))
# 
# ChlWinterDF$z[is.na(ChlWinterDF$z)] = 0
# newWin = regrid(ChlWinterDF,n1=length(newx),n2=length(newy),method="idw")
# ChlSpringDF$z[is.na(ChlSpringDF$z)] = 0
# newSpr = regrid(ChlSpringDF,n1=length(newx),n2=length(newy),method="idw")
# ChlSummerDF$z[is.na(ChlSummerDF$z)] = 0
# newSum = regrid(ChlSummerDF,n1=length(newx),n2=length(newy),method="idw")
# ChlFallDF$z[is.na(ChlFallDF$z)] = 0
# newFall = regrid(ChlFallDF,n1=length(newx),n2=length(newy),method="idw")
# 
# predictWinter$Chl0 = unlist(stack(data.frame(matrix(newWin$z,ncol=length(newx),byrow=TRUE)))[,1])
# predictSpring$Chl0 = unlist(stack(data.frame(matrix(newSpr$z,ncol=length(newx),byrow=TRUE)))[,1])
# predictSummer$Chl0 = unlist(stack(data.frame(matrix(newSum$z,ncol=length(newx),byrow=TRUE)))[,1])
# predictFall$Chl0 = unlist(stack(data.frame(matrix(newFall$z,ncol=length(newx),byrow=TRUE)))[,1])

# save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/SeasonalPredictionData.Rdata',sep=""))
save(predict_1,predict_2,predict_3,predict_4,predict_5,predict_6,
     predict_7,predict_8,predict_9,predict_10,predict_11,predict_12,
     file=paste(outDir,'/MonthlyPredictionData.Rdata',sep=""))

# Eddy Dist -----------------------------------------------------

inDir = 'E:/Chpt_3/Eddies/Grids'

# anticyclonic eddies
AfileList = dir(inDir,"AEddy",full.names=TRUE,recursive=FALSE)
# AEddy_Winter = numeric()
# AEddy_Spring = numeric()
# AEddy_Summer = numeric()
# AEddy_Fall = numeric()

for (k in 1:length(AfileList)){
  # get 6-digit datestamps from file names
  fileDate = str_extract(AfileList[k],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  thisTime = as_datetime(fileDate,format='%Y%m%d',tz="UTC")
  thisMonth = month(thisTime)
  
  # if (month(thisTime)==1){
  #   load(AfileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   AEddy_Winter = rbind(AEddy_Winter,dataVec)
  # } else if (month(thisTime)==4){
  #   load(AfileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   AEddy_Spring = rbind(AEddy_Spring,dataVec)  
  # } else if (month(thisTime)==7){
  #   load(AfileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   AEddy_Summer = rbind(AEddy_Summer,dataVec)
  # } else if (month(thisTime)==10){
  #   load(AfileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   AEddy_Fall = rbind(AEddy_Fall,dataVec)
  # }
  
  dataVec = stack(data.frame(data))[,1]
  if (!exists(eval(parse(text=paste('\'AEddy_',thisMonth,'\'',sep=""))))){
    eval(parse(text=paste('AEddy_',thisMonth,' = dataVec',sep="")))
  }else{
    eval(parse(text=paste('AEddy_',thisMonth,' = rbind(AEddy_',thisMonth,',dataVec)',sep="")))
  }
  
}

# predictWinter$AEddyDist0 = apply(AEddy_Winter,2,mean,na.rm=TRUE)
# predictSpring$AEddyDist0 = apply(AEddy_Spring,2,mean,na.rm=TRUE)
# predictSummer$AEddyDist0 = apply(AEddy_Summer,2,mean,na.rm=TRUE)
# predictFall$AEddyDist0 = apply(AEddy_Fall,2,mean,na.rm=TRUE)

for (l in 1:12){
  eval(parse(text=paste('predict_',l,' = data.frame(AEddyDist0 = apply(AEddy_',l,',MARGIN=2,mean,na.rm=TRUE))',sep="")))
}

# Cyclonic eddies
CfileList = dir(inDir,"CEddy",full.names=TRUE,recursive=FALSE)

# CEddy_Winter = numeric()
# CEddy_Spring = numeric()
# CEddy_Summer = numeric()
# CEddy_Fall = numeric()

for (k in 1:length(CfileList)){
  # get 6-digit datestamps from file names
  fileDate = str_extract(CfileList[k],"\\d\\d\\d\\d\\d\\d\\d\\d") 
  thisTime = as_datetime(fileDate,format='%Y%m%d',tz="UTC")
  thisMonth = month(thisTime)
  
  # if (month(thisTime)==1){
  #   load(CfileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   CEddy_Winter = rbind(CEddy_Winter,dataVec)
  # } else if (month(thisTime)==4){
  #   load(CfileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   CEddy_Spring = rbind(CEddy_Spring,dataVec)  
  # } else if (month(thisTime)==7){
  #   load(CfileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   CEddy_Summer = rbind(CEddy_Summer,dataVec)
  # } else if (month(thisTime)==10){
  #   load(CfileList[k])
  #   dataVec = stack(data.frame(data))[,1]
  #   CEddy_Fall = rbind(CEddy_Fall,dataVec)
  # }
  
  dataVec = stack(data.frame(data))[,1]
  if (!exists(eval(parse(text=paste('\'CEddy_',thisMonth,'\'',sep=""))))){
    eval(parse(text=paste('CEddy_',thisMonth,' = dataVec',sep="")))
  }else{
    eval(parse(text=paste('CEddy_',thisMonth,' = rbind(CEddy_',thisMonth,',dataVec)',sep="")))
  }
  
}

# predictWinter$CEddyDist0 = apply(CEddy_Winter,2,mean,na.rm=TRUE)
# predictSpring$CEddyDist0 = apply(CEddy_Spring,2,mean,na.rm=TRUE)
# predictSummer$CEddyDist0 = apply(CEddy_Summer,2,mean,na.rm=TRUE)
# predictFall$CEddyDist0 = apply(CEddy_Fall,2,mean,na.rm=TRUE)

for (l in 1:12){
  eval(parse(text=paste('predict_',l,' = data.frame(CEddyDist0 = apply(CEddy_',l,',MARGIN=2,mean,na.rm=TRUE))',sep="")))
}

# save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/SeasonalPredictionData.Rdata',sep=""))
save(predict_1,predict_2,predict_3,predict_4,predict_5,predict_6,
     predict_7,predict_8,predict_9,predict_10,predict_11,predict_12,
     file=paste(outDir,'/MonthlyPredictionData.Rdata',sep=""))

# Convert prediction data to data frames -----------------

# predictWinter = data.frame(predictWinter)
# predictSpring = data.frame(predictSpring)
# predictSummer = data.frame(predictSummer)
# predictFall = data.frame(predictFall)

# rename columns to match model covars
# colnames(predictWinter) = c("lat","lon","EKE0","SSH0","Sal0","Sal200","Sal400","Sal700",
#                             "Temp0","Temp200","Temp400","Temp700","VelAsp0","VelAsp200",
#                             "VelAsp400","VelAsp700","VelMag0","VelMag200","VelMag400",
#                             "VelMag700","FSLE0","Chl0","AEddyDist0","CEddyDist0")
# colnames(predictSpring) = c("lat","lon","EKE0","SSH0","Sal0","Sal200","Sal400","Sal700",
#                             "Temp0","Temp200","Temp400","Temp700","VelAsp0","VelAsp200",
#                             "VelAsp400","VelAsp700","VelMag0","VelMag200","VelMag400",
#                             "VelMag700","FSLE0","Chl0","AEddyDist0","CEddyDist0")
# colnames(predictSummer) = c("lat","lon","EKE0","SSH0","Sal0","Sal200","Sal400","Sal700",
#                             "Temp0","Temp200","Temp400","Temp700","VelAsp0","VelAsp200",
#                             "VelAsp400","VelAsp700","VelMag0","VelMag200","VelMag400",
#                             "VelMag700","FSLE0","Chl0","AEddyDist0","CEddyDist0")
# colnames(predictFall) = c("lat","lon","EKE0","SSH0","Sal0","Sal200","Sal400","Sal700",
#                             "Temp0","Temp200","Temp400","Temp700","VelAsp0","VelAsp200",
#                             "VelAsp400","VelAsp700","VelMag0","VelMag200","VelMag400",
#                             "VelMag700","FSLE0","Chl0","AEddyDist0","CEddyDist0")
# colnames(predictWinter) = c("lat","lon","EKE0","SSH0","Sal0","Temp0","VelAsp0","VelMag0","FSLE0","Chl0","AEddyDist0","CEddyDist0")
# colnames(predictSpring) = c("lat","lon","EKE0","SSH0","Sal0","Temp0","VelAsp0","VelMag0","FSLE0","Chl0","AEddyDist0","CEddyDist0")
# colnames(predictSummer) = c("lat","lon","EKE0","SSH0","Sal0","Temp0","VelAsp0","VelMag0","FSLE0","Chl0","AEddyDist0","CEddyDist0")
# colnames(predictFall) = c("lat","lon","EKE0","SSH0","Sal0","Temp0","VelAsp0","VelMag0","FSLE0","Chl0","AEddyDist0","CEddyDist0")

# transform data as necessary
# # -Inf will occur in Chl0 + FSLE0 due to 0 values, set to NA
# predictWinter$Chl0[predictWinter$Chl0>10] = 10 # cap crazy high values near coast
# predictWinter$log_Chl0 = log10(predictWinter$Chl0)
# predictWinter$log_Chl0[is.infinite(predictWinter$log_Chl0)] = NA
# predictWinter$log_abs_FSLE0 = log10(abs(predictWinter$FSLE0))
# predictWinter$log_abs_FSLE0[is.infinite(predictWinter$log_abs_FSLE0)] = NA
# predictWinter$sqrt_CEddyDist0 = sqrt(predictWinter$CEddyDist0)
# predictWinter$sqrt_AEddyDist0 = sqrt(predictWinter$AEddyDist0)
# predictWinter$sqrt_VelAsp0 = sqrt(predictWinter$VelAsp0)
# predictWinter$sqrt_VelAsp700 = sqrt(predictWinter$VelAsp700)
# predictWinter$sqrt_EKE0 = sqrt(predictWinter$EKE0)
# 
# predictSpring$Chl0[predictSpring$Chl0>10] = 10
# predictSpring$log_Chl0 = log10(predictSpring$Chl0)
# predictSpring$log_Chl0[is.infinite(predictSpring$log_Chl0)] = NA
# predictSpring$log_abs_FSLE0 = log10(abs(predictSpring$FSLE0))
# predictSpring$log_abs_FSLE0[is.infinite(predictSpring$log_abs_FSLE0)] = NA
# predictSpring$sqrt_CEddyDist0 = sqrt(predictSpring$CEddyDist0)
# predictSpring$sqrt_AEddyDist0 = sqrt(predictSpring$AEddyDist0)
# predictSpring$sqrt_VelAsp0 = sqrt(predictSpring$VelAsp0)
# predictSpring$sqrt_VelAsp700 = sqrt(predictSpring$VelAsp700)
# predictSpring$sqrt_EKE0 = sqrt(predictSpring$EKE0)
# 
# predictSummer$Chl0[predictSummer$Chl0>10] = 10
# predictSummer$log_Chl0 = log10(predictSummer$Chl0)
# predictSummer$log_Chl0[is.infinite(predictSummer$log_Chl0)] = NA
# predictSummer$log_abs_FSLE0 = log10(abs(predictSummer$FSLE0))
# predictSummer$log_abs_FSLE0[is.infinite(predictSummer$log_abs_FSLE0)] = NA
# predictSummer$sqrt_CEddyDist0 = sqrt(predictSummer$CEddyDist0)
# predictSummer$sqrt_AEddyDist0 = sqrt(predictSummer$AEddyDist0)
# predictSummer$sqrt_VelAsp0 = sqrt(predictSummer$VelAsp0)
# predictSummer$sqrt_VelAsp700 = sqrt(predictSummer$VelAsp700)
# predictSummer$sqrt_EKE0 = sqrt(predictSummer$EKE0)
# 
# predictFall$Chl0[predictFall$Chl0>10] = 10
# predictFall$log_Chl0 = log10(predictFall$Chl0)
# predictFall$log_Chl0[is.infinite(predictFall$log_Chl0)] = NA
# predictFall$log_abs_FSLE0 = log10(abs(predictFall$FSLE0))
# predictFall$log_abs_FSLE0[is.infinite(predictFall$log_abs_FSLE0)] = NA
# predictFall$sqrt_CEddyDist0 = sqrt(predictFall$CEddyDist0)
# predictFall$sqrt_AEddyDist0 = sqrt(predictFall$AEddyDist0)
# predictFall$sqrt_VelAsp0 = sqrt(predictFall$VelAsp0)
# predictFall$sqrt_VelAsp700 = sqrt(predictFall$VelAsp700)
# predictFall$sqrt_EKE0 = sqrt(predictFall$EKE0)

# get elevation values to later remove prediction values on land
load('E:/Chpt_3/GEBCO/DownsampledGrid_08deg.Rdata')
x = seq(278,297,by=0.08)
y = seq(24,44,by=0.08)
depth = unlist(stack(data.frame(matrix(newGrid$Depth,ncol=length(x),byrow=TRUE)))[,1])
depth[depth>=0] = NA
# predictWinter$Depth = depth
# predictSpring$Depth = depth
# predictSummer$Depth = depth
# predictFall$Depth = depth

# # set margin values = NA
# lim300 = which(predictWinter$Depth>-300)
# var200 = which(str_detect(colnames(predictWinter),"200"))
# predictWinter[lim300,var200] = NA
# predictSpring[lim300,var200] = NA
# predictSummer[lim300,var200] = NA
# predictFall[lim300,var200] = NA
# 
# lim500 = which(predictWinter$Depth>-500)
# var400 = which(str_detect(colnames(predictWinter),"400"))
# predictWinter[lim500,var400] = NA
# predictSpring[lim500,var400] = NA
# predictSummer[lim500,var400] = NA
# predictFall[lim500,var400] = NA
# 
# lim800 = which(predictWinter$Depth>-800)
# var700 = which(str_detect(colnames(predictWinter),"700"))
# predictWinter[lim800,var700] = NA
# predictSpring[lim800,var700] = NA
# predictSummer[lim800,var700] = NA
# predictFall[lim800,var700] = NA

cnames = c("lat","lon","EKE0","SSH0","Sal0","Temp0","VelAsp0","VelMag0","FSLE0","Chl0","AEddyDist0","CEddyDist0")

for (i in 1:12){
  predict_l = data.frame(predict_l)
  colnames(predict_l) = cnames
  eval(parse(text=paste('predict_',l,'$Depth = depth',sep="")))
}

# save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))
save(predict_1,predict_2,predict_3,predict_4,predict_5,predict_6,
     predict_7,predict_8,predict_9,predict_10,predict_11,predict_12,
     file=paste(outDir,'/MonthlyPredictionData.Rdata',sep=""))
