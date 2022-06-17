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
      } else if (month(thisTime)==4){
        load(fileList[varFiles[k]])
        dataVec = stack(data.frame(data))[,1]
        eval(parse(text=paste('thisVar_Spring = rbind(thisVar_Spring,dataVec)',sep="")))
      } else if (month(thisTime)==7){
        load(fileList[varFiles[k]])
        dataVec = stack(data.frame(data))[,1]
        eval(parse(text=paste('thisVar_Summer = rbind(thisVar_Summer,dataVec)',sep="")))
      } else if (month(thisTime)==10){
        load(fileList[varFiles[k]])
        dataVec = stack(data.frame(data))[,1]
        eval(parse(text=paste('thisVar_Fall = rbind(thisVar_Fall,dataVec)',sep="")))
      }
      
    }
    
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

for (k in 1:length(fileList)){
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
bChlSpringDF = data.frame(z=stack(data.frame(apply(Chl_Spring,2,mean,na.rm=TRUE)))[,1],
                         y=rep(lats,length.out=length(lons)*length(lats)),
                         x=rep(lons,each=length(lats)))
ChlSummerDF = data.frame(z=stack(data.frame(apply(Chl_Summer,2,mean,na.rm=TRUE)))[,1],
                         y=rep(lats,length.out=length(lons)*length(lats)),
                         x=rep(lons,each=length(lats)))
ChlFallDF = data.frame(z=stack(data.frame(apply(Chl_Fall,2,mean,na.rm=TRUE)))[,1],
                         y=rep(lats,length.out=length(lons)*length(lats)),
                         x=rep(lons,each=length(lats)))

ChlWinterDF$z[is.na(ChlWinterDF$z)] = -1
newWin = regrid(ChlWinterDF,n1=length(newx),n2=length(newy),method="idw")
ChlSpringDF$z[is.na(ChlSpringDF$z)] = -1
newSpr = regrid(ChlSpringDF,n1=length(newx),n2=length(newy),method="idw")
ChlSummerDF$z[is.na(ChlSummerDF$z)] = -1
newSum = regrid(ChlSummerDF,n1=length(newx),n2=length(newy),method="idw")
ChlFallDF$z[is.na(ChlFallDF$z)] = -1
newFall = regrid(ChlFallDF,n1=length(newx),n2=length(newy),method="idw")

predictWinter$Chl0 = unlist(newWin$z)
predictSpring$Chl0 = unlist(newSpr$z)
predictSummer$Chl0 = unlist(newSum$z)
predictFall$Chl0 = unlist(newFall$z)

save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))

# Eddy Dist -----------------------------------------------------


save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))

# Convert prediction data to data frames -----------------

predictWinter = data.frame(predictWinter)
predictSpring = data.frame(predictSpring)
predictSummer = data.frame(predictSummer)
predictFall = data.frame(predictFall)

# rename columns to match model covars

# get continental mask and make all data points on land = NA


save(predictWinter,predictSummer,predictSpring,predictFall,file=paste(outDir,'/PredictionData.Rdata',sep=""))
