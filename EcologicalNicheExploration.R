
inDir = 'J:/Chpt_3/ModelData'
fileList = list.files(path=inDir,pattern=paste('_masterDF.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)

for (i in 1:length(fileList)){
  
  spec = str_remove(str_remove(fileList[i],'_masterDF.csv'),paste(inDir,'/',sep=""))
  data = data.frame(read.csv(fileList[i]))
  
  # TO DO: find bins with presence (don't plot all the 0's)
  
  optTemp = sum(data$Pres*data$Temp0)/sum(data$Pres)
  
  plot(data$Temp0,data$Pres,type="p",main=spec)
  abline(v=optTemp,col="red",lty=2)
  
  optDist = data$Temp0-optTemp
  
  plot(abs(optDist),data$Pres,type="p",main=spec)
  #Add a fit line to this plot
}