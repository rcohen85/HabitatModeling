# Load regional models and use seasonal climatology covar grids to predict seasonal presence

library(stringr)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(viridis)
library(maps)

modDir = 'J:/Chpt_3/GAM_Output'
predDir = 'E:/Chpt_3/Predictions'
lon = seq(278,297,by=0.08)-360
lat = seq(24,44,by=0.08)

species = list.dirs(modDir,recursive=FALSE)
load(paste(predDir,'/PredictionData.Rdata',sep=""))
season = c("Fall","Winter","Spring","Summer")

# Plot prediction data -----------------------------
# for (i in 3:dim(predictWinter)[2]){
# 
  # thisVar = names(predictFall)[i]
  # varMin = min(c(predictFall[[i]],predictWinter[[i]],predictSpring[[i]],predictSummer[[i]]),na.rm=TRUE)
  # varMax = max(c(predictFall[[i]],predictWinter[[i]],predictSpring[[i]],predictSummer[[i]]),na.rm=TRUE)
  # 
  # 
  # Fall = ggplot(predictFall,aes(x=lon-360,y=lat)
  #        )+geom_tile(aes(fill=predictFall[[i]])
  #        )+scale_fill_viridis(limits=c(varMin,varMax)
  #        )+guides(fill = guide_colorbar(title="Fall")
  #        )+geom_point(data=HARPs,aes(x=lon, y=lat),color="#FDFEFE",size=2,shape=19)
  # Winter = ggplot(predictWinter,aes(x=lon-360,y=lat)
  #       )+geom_tile(aes(fill=predictWinter[[i]])
  #       )+scale_fill_viridis(limits=c(varMin,varMax)
  #       )+guides(fill = guide_colorbar(title="Winter")
  #       )+geom_point(data=HARPs,aes(x=lon, y=lat),color="#FDFEFE",size=2,shape=19)
  # Spring = ggplot(predictSpring,aes(x=lon-360,y=lat)
  #       )+geom_tile(aes(fill=predictSpring[[i]])
  #       )+scale_fill_viridis(limits=c(varMin,varMax)
  #       )+guides(fill = guide_colorbar(title="Spring")
  #       )+geom_point(data=HARPs,aes(x=lon, y=lat),color="#FDFEFE",size=2,shape=19)
  # Summer = ggplot(predictSummer,aes(x=lon-360,y=lat)
  #       )+geom_tile(aes(fill=predictSummer[[i]])
  #       )+scale_fill_viridis(limits=c(varMin,varMax)
  #       )+guides(fill = guide_colorbar(title="Summer")
  #       )+geom_point(data=HARPs,aes(x=lon, y=lat),color="#FDFEFE",size=2,shape=19)
# 
#   png(file=paste(predDir,"/PredictionData/",thisVar,".png",sep=""),height=700,width=700,units="px",res=100)
  # grid.arrange(Winter,Spring,Summer,Fall,ncol=2,nrow=2,top=thisVar,layout_matrix=rbind(c(1,2),c(4,3)))
#   while (dev.cur()>1) {dev.off()}
# 
#   # ggsave(filename=paste(predDir,"/PredictionData/",season[k],"_",thisVar,".png",sep=""),
#   #        device="png",
#   #        width=300,
#   #        height=300,
#   #        units="px",
#   #        scale=7)
# 
# }

# DEdelta = matrix(nrow=length(species),ncol=2)

# Predict and plot ---------------------------

HARPs = data.frame(t(data.frame(c(41.06165, -66.35155), # WAT_HZ
                                c(40.26333,-67.9861),  # WAT_OC
                                c(39.83295, -69.98194),  # WAT_NC
                                c(39.19192, -72.22735),  # WAT_BC
                                c(38.37337, -73.36985),  # WAT_WC
                                c(37.16452, -74.46585),  # NFC
                                c(35.30183,-74.8789),  # HAT_A & HAT_B
                                c(33.66992, -75.9977),   # WAT_GS
                                c(32.10527, -77.09067),  # WAT_BP
                                c(30.58295, -77.39002))))  # WAT_BS
colnames(HARPs) = c("lat","lon")
rownames(HARPs) = c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS")

# Center & scale prediction data
predictWinter[,3:12] = scale(predictWinter[,3:12],center=TRUE,scale=TRUE)
predictSpring[,3:12] = scale(predictSpring[,3:12],center=TRUE,scale=TRUE)
predictSummer[,3:12] = scale(predictSummer[,3:12],center=TRUE,scale=TRUE)
predictFall[,3:12] = scale(predictFall[,3:12],center=TRUE,scale=TRUE)

for (i in 1:length(species)){
  
  specName = str_remove(species[i],paste(modDir,'/',sep=""))
  load(paste(species[i],"/",specName,'_WeeklyRegionalModel.Rdata',sep=""))
  load(paste(species[i],"/",specName,'_AlternateModel.Rdata',sep=""))
  
  DEOpt = ((optWeekMod$null.deviance-optWeekMod$deviance)/optWeekMod$null.deviance)*100
  # DEAlt = ((altMod$null.deviance-altMod$deviance)/altMod$null.deviance)*100
  # DEchange = DEOpt - DEAlt
  # 
  # DEdelta[i,] = cbind(specName,round(DEchange,digits=1))
  
  # predict presence
  winPred = predict.gam(topMods[['96']],predictWinter,type="response",na.action=na.pass)
  sprPred = predict.gam(topMods[['96']],predictSpring,type="response",na.action=na.pass)
  sumPred = predict.gam(topMods[['96']],predictSummer,type="response",na.action=na.pass)
  fallPred = predict.gam(topMods[['96']],predictFall,type="response",na.action=na.pass)
  
  winPred = predict(optWeekMod,predictWinter,full=TRUE,type="response",backtransform=FALSE)
  sprPred = predict(optWeekMod,predictSpring,full=TRUE,type="response",backtransform=FALSE)
  sumPred = predict(optWeekMod,predictSummer,full=TRUE,type="response",backtransform=FALSE)
  fallPred = predict(optWeekMod,predictFall,full=TRUE,type="response",backtransform=FALSE)
  
  winPred = predict(optWeekMod,predictWinter,full=TRUE,type="link",backtransform=TRUE)
  sprPred = predict(optWeekMod,predictSpring,full=TRUE,type="link",backtransform=TRUE)
  sumPred = predict(optWeekMod,predictSummer,full=TRUE,type="link",backtransform=TRUE)
  fallPred = predict(optWeekMod,predictFall,full=TRUE,type="link",backtransform=TRUE)
  
  # winPred = predict.gam(optWeekMod,predictWinter,type="response",na.action=na.pass)
  # sprPred = predict.gam(optWeekMod,predictSpring,type="response",na.action=na.pass)
  # sumPred = predict.gam(optWeekMod,predictSummer,type="response",na.action=na.pass)
  # fallPred = predict.gam(optWeekMod,predictFall,type="response",na.action=na.pass)
  
  # # undo link function on predictions, if predictions are type="link"
  # winPred = exp(winPred)
  # sprPred = exp(sprPred)
  # sumPred = exp(sumPred)
  # fallPred = exp(fallPred)
  
  
  if (specName=="Blainville"){ # clip presence to max actually possible per week for Md
  # winPred[winPred>2016] = 2016
  # sprPred[sprPred>2016] = 2016
  # sumPred[sumPred>2016] = 2016
  # fallPred[fallPred>2016] = 2016
  } else if (specName=="Kogia"){ # cap these species to improve dynamic range of color scale
    winPred[winPred>30] = 30
    sprPred[sprPred>30] = 30
    sumPred[sumPred>30] = 30
    fallPred[fallPred>30] = 30
  } else if (specName=="True"){ 
    winPred[winPred>20] = 20
    sprPred[sprPred>20] = 20
    sumPred[sumPred>20] = 20
    fallPred[fallPred>20] = 20
  } else if (specName=="SBCD"){ 
    winPred[winPred>1000] = 1000
    sprPred[sprPred>1000] = 1000
    sumPred[sumPred>1000] = 1000
    fallPred[fallPred>1000] = 1000
  } else if (specName=="Risso"){
    winPred[winPred>450] = 450
    sprPred[sprPred>450] = 450
    sumPred[sumPred>450] = 450
    fallPred[fallPred>450] = 450
  } else if (specName=="Cuvier"){ 
    winPred[winPred>350] = 350
    sprPred[sprPred>350] = 350
    sumPred[sumPred>350] = 350
    fallPred[fallPred>350] = 350
  } else if (specName=="SpermWhale"){ 
    winPred[winPred>800] = 800
    sprPred[sprPred>800] = 800
    sumPred[sumPred>800] = 800
    fallPred[fallPred>800] = 800
  }
  
  # create continental mask
  mask = rep(1,dim(predictWinter)[1])
  mask[is.na(predictWinter$Depth)] = 0
  if (specName!="SBCD" & specName!="Risso"){ # get rid of shallow predictions for deep divers
  mask[predictWinter$Depth>-200] = 0}

  # set land = NA
  winPred[which(mask==0)] = NA
  sprPred[which(mask==0)] = NA
  sumPred[which(mask==0)] = NA
  fallPred[which(mask==0)] = NA
  
  # # get rid of any predictions in cells not deep enough for the covars in this model
  # thisForm = as.character(optWeekMod$formula)[3]
  # startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
  # termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
  # termInd = c(0,termInd,str_length(thisForm)+1)
  # allTerms = character()
  # for (j in 1:length(termInd)-1){
  #   thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
  #   allTerms = c(allTerms,thisTerm)
  # }
  # nullTerm = str_which(allTerms,"1")
  # allTerms = allTerms[-nullTerm]
  # if (any(str_detect(allTerms,"700"))){
  #   dLim = -800
  # } else if (any(str_detect(allTerms,"400"))){
  #   dLim = -500
  # } else if (any(str_detect(allTerms,"200"))){
  #   dLim = -300
  # }
# 
#   winPred[which(predictWinter$Depth>dLim)] = NA
#   sprPred[which(predictSpring$Depth>dLim)] = NA
#   sumPred[which(predictSummer$Depth>dLim)] = NA
#   fallPred[which(predictFall$Depth>dLim)] = NA

  # make consistent fill limits
  low = min(c(winPred,sprPred,sumPred,fallPred),na.rm=TRUE)
  # high = min(c(2016,max(c(winPred,sprPred,sumPred,fallPred),na.rm=TRUE)))
  high = max(c(winPred,sprPred,sumPred,fallPred),na.rm=TRUE)
  
  # make data frames for plotting
  WinDF = data.frame(lat=rep(lat,length.out=length(lon)*length(lat)),
                     lon=rep(lon,each=length(lat)),z=winPred)
  SprDF = data.frame(lat=rep(lat,length.out=length(lon)*length(lat)),
                     lon=rep(lon,each=length(lat)),z=sprPred)
  SumDF = data.frame(lat=rep(lat,length.out=length(lon)*length(lat)),
                     lon=rep(lon,each=length(lat)),z=sumPred)
  FallDF = data.frame(lat=rep(lat,length.out=length(lon)*length(lat)),
                      lon=rep(lon,each=length(lat)),z=fallPred)
  
  # create plots
  winPlot = ggplot(WinDF,aes(x=lon,y=lat)
  )+geom_tile(aes(fill=z)
  )+geom_point(data=HARPs,aes(x=lon, y=lat),color="#FDFEFE",size=2,shape=19
  )+scale_fill_viridis(limits=c(low,high)
  )+ggtitle('Winter'
  )+theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_blank())
  
  sprPlot = ggplot(SprDF,aes(x=lon,y=lat)
  )+geom_tile(aes(fill=z)
  )+geom_point(data=HARPs,aes(x=lon, y=lat),color="#FDFEFE",size=2,shape=19
  )+scale_fill_viridis(limits=c(low,high)
  )+ggtitle('Spring'
  )+theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_blank())
  
  sumPlot = ggplot(SumDF,aes(x=lon,y=lat)
  )+geom_tile(aes(fill=z)
  )+geom_point(data=HARPs,aes(x=lon, y=lat),color="#FDFEFE",size=2,shape=19
  )+scale_fill_viridis(limits=c(low,high)
  )+ggtitle('Summer'
  )+theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_blank())
  
  fallPlot = ggplot(FallDF,aes(x=lon,y=lat)
  )+geom_tile(aes(fill=z)
  )+geom_point(data=HARPs,aes(x=lon, y=lat),color="#FDFEFE",size=2,shape=19
  )+scale_fill_viridis(limits=c(low,high)
  )+guides(
  )+ggtitle('Fall'
  )+theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_blank())
  
  # if (i==1){ # plot coastline on its own to later overlay on top of prediction maps
  #   Coast = ggplot(FallDF,aes(x=lon,y=lat)
  #   )+borders(xlim=c(min(lon),max(lon)),
  #             ylim=c(min(lat),max(lat)),
  #             fill=NA,colour="black"
  #   )+coord_cartesian(xlim=c(min(lon),max(lon)),
  #                     ylim=c(min(lat),max(lat)),
  #   )+theme(panel.grid.major.x=element_blank(),
  #           panel.grid.minor.x=element_blank()
  #   )+theme_minimal()
  #   ggsave(Coast,file=paste(predDir,'/',"Coastline.pdf",sep=""),device="pdf")
  # }
  
  # arrange plots in grid and save
  png(file=paste(predDir,'/',specName,'.png',sep=""),width=750,height=800)
  grid.arrange(winPlot,sprPlot,fallPlot,sumPlot,nrow=2,ncol=2,top=specName)
  while (dev.cur()>1) {dev.off()}
  
  # pdf(file=paste(predDir,'/',specName,'.pdf',sep=""))
  # grid.arrange(winPlot,sprPlot,fallPlot,sumPlot,nrow=2,ncol=2,top=specName)
  # while (dev.cur()>1) {dev.off()}
  
  
}
