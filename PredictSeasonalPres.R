# Load regional models and use seasonal climatology covar grids to predict seasonal presence

library(stringr)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(viridis)

modDir = 'J:/Chpt_3/GAM_Output'
predDir = 'J:/Chpt_3/Predictions'
lon = seq(278,297,by=0.08)
lat = seq(24,44,by=0.08)

species = list.dirs(modDir,recursive=FALSE)
load(paste(predDir,'/PredictionData.Rdata',sep=""))
datList = list(predictFall,predictWinter,predictSpring,predictSummer)
season = c("Fall","Winter","Spring","Summer")

# for (k in 1:length(datList)){
for (i in 3:dim(datList[[k]])[2]){
  
  thisVar = names(predictFall)[i]
  varMin = min(c(predictFall[[i]],predictWinter[[i]],predictSpring[[i]],predictSummer[[i]]),na.rm=TRUE)
  varMax = max(c(predictFall[[i]],predictWinter[[i]],predictSpring[[i]],predictSummer[[i]]),na.rm=TRUE)
  
  
  Fall = ggplot(predictFall,aes(x=lon,y=lat)
         )+geom_tile(aes(fill=predictFall[[i]])
         )+scale_fill_viridis(limits=c(varMin,varMax)
         )+guides(fill = guide_colorbar(title="Fall"))
  Winter = ggplot(predictWinter,aes(x=lon,y=lat)
        )+geom_tile(aes(fill=predictWinter[[i]])
        )+scale_fill_viridis(limits=c(varMin,varMax)
        )+guides(fill = guide_colorbar(title="Winter"))
  Spring = ggplot(predictSpring,aes(x=lon,y=lat)
        )+geom_tile(aes(fill=predictSpring[[i]])
        )+scale_fill_viridis(limits=c(varMin,varMax)
        )+guides(fill = guide_colorbar(title="Spring"))
  Summer = ggplot(predictSummer,aes(x=lon,y=lat)
        )+geom_tile(aes(fill=predictSummer[[i]])
        )+scale_fill_viridis(limits=c(varMin,varMax)
        )+guides(fill = guide_colorbar(title="Summer"))
  
  png(file=paste(predDir,"/PredictionData/",thisVar,"Clipped.png",sep=""),height=700,width=700,units="px",res=100)
  grid.arrange(Fall,Winter,Spring,Summer,ncol=2,nrow=2,top=thisVar,layout_matrix=rbind(c(1,2),c(4,3)))
  while (dev.cur()>1) {dev.off()}
  
  # ggsave(filename=paste(predDir,"/PredictionData/",season[k],"_",thisVar,".png",sep=""),
  #        device="png",
  #        width=300,
  #        height=300,
  #        units="px",
  #        scale=7)
  
}
# }

for (i in 1:length(species)){
  
  load(paste(species[i],'/WeeklyRegionalModel.Rdata',sep=""))
  
  specName = str_remove(species[i],paste(modDir,'/',sep=""))
  
  # create continental mask
  mask = rep(1,dim(predictWinter)[1])
  mask[is.na(predictWinter$Depth)] = 0
  
  # predict presence
  winPred = predict.gam(optWeekMod,predictWinter,type="response")
  sprPred = predict.gam(optWeekMod,predictSpring,type="response")
  sumPred = predict.gam(optWeekMod,predictSummer,type="response")
  fallPred = predict.gam(optWeekMod,predictFall,type="response")

  # set land = NA
  winPred[which(mask==0)] = NA
  sprPred[which(mask==0)] = NA
  sumPred[which(mask==0)] = NA
  fallPred[which(mask==0)] = NA
  
  # make consistent fill limits
  low = min(c(winPred,sprPred,sumPred,fallPred),na.rm=TRUE)
  high = min(c(2016,max(c(winPred,sprPred,sumPred,fallPred),na.rm=TRUE)))
  
  # make data frames for plotting
  WinDF = data.frame(y=rep(lat,length.out=length(lon)*length(lat)),
                     x=rep(lon,each=length(lat)),z=winPred)
  SprDF = data.frame(y=rep(lat,length.out=length(lon)*length(lat)),
                     x=rep(lon,each=length(lat)),z=sprPred)
  SumDF = data.frame(y=rep(lat,length.out=length(lon)*length(lat)),
                     x=rep(lon,each=length(lat)),z=sumPred)
  FallDF = data.frame(y=rep(lat,length.out=length(lon)*length(lat)),
                      x=rep(lon,each=length(lat)),z=fallPred)
  
  # create plots
  winPlot = ggplot(WinDF,aes(x=x,y=y)
  )+geom_tile(aes(fill=z)
  )+scale_fill_viridis(limits=c(low,high)
  )+ggtitle('Winter')
  
  sprPlot = ggplot(SprDF,aes(x=x,y=y)
  )+geom_tile(aes(fill=z)
  )+scale_fill_viridis(limits=c(low,high)
  )+ggtitle('Spring')
  
  sumPlot = ggplot(SumDF,aes(x=x,y=y)
  )+geom_tile(aes(fill=z)
  )+scale_fill_viridis(limits=c(low,high)
  )+ggtitle('Summer')
  
  fallPlot = ggplot(FallDF,aes(x=x,y=y)
  )+geom_tile(aes(fill=z)
  )+scale_fill_viridis(limits=c(low,high)
  )+ggtitle('Fall')
  
  # arrange plots in grid and save
  png(file=paste(predDir,'/',specName,'.png',sep=""),width=750,height=800,)
  grid.arrange(winPlot,sprPlot,fallPlot,sumPlot,nrow=2,ncol=2,top=specName)
  while (dev.cur()>1) {dev.off()}
  
  
}