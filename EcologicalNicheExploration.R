library(stringr)
inDir = 'E:/ModelingCovarData/Master_DFs'
outDir = 'E:/ModelingCovarData/Ecological Niche Plots/Sowerby'
fileList = list.files(path=inDir,pattern=paste('_masterDF.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)
# covarList = c("Temp0","Temp100","Temp200","Temp300","Temp400","Temp500","Temp600","Temp700")
# covarList = c("Sal0","Sal100","Sal200","Sal300","Sal400","Sal500","Sal600","Sal700")
# covarList = c("SSH0")
# covarList = c("Chl0")
covarList = c("FSLE0")

for (i in 1:length(fileList)){
  
  spec = str_remove(str_remove(fileList[i],'_masterDF.csv'),paste(inDir,'/',sep=""))
  data = data.frame(read.csv(fileList[i]))
  
  # TO DO: find bins with presence (don't plot all the 0's)
  data$Pres = round(data$Pres)
  data = subset(data, data$Pres>0)
  # data = subset(data, data$Site=="WC")
  
  # attach(mtcars)
  # par(mfrow=c(8,2))
  
  for (j in 1:length(covarList)) {
  
  # jpeg(paste(outDir,"/",spec,"_",covarList[j],"_allSites.jpg",sep=""))
  jpeg(paste(outDir,"/",spec,"_",covarList[j],"_",data$Site[1],".jpg",sep=""))
    
  attach(mtcars)
  par(mfrow=c(1,2))
  optCovar = sum(data$Pres*data[[covarList[j]]])/sum(data$Pres)
  
  # plot(data[[covarList[j]]],data$Pres,type="p",main=paste(spec,'&',covarList[j]))
  plot(data[[covarList[j]]],data$Pres,type="p",main=paste(spec,'at',data$Site[1]),
       xlab=paste(covarList[j]))
  abline(v=optCovar,col="red",lty=2)
  
  optDist = data[[covarList[j]]]-optCovar
  
  # plot(abs(optDist),data$Pres,type="p",main=paste(spec,'&',covarList[j]))
  plot(abs(optDist),data$Pres,type="p",main=paste(spec,'at',data$Site[1]))
  #Add a fit line to this plot
  abline(lm(data$Pres~data[[covarList[j]]]),col="red")
  
  dev.off()

  }

}
