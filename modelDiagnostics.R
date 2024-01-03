# Calculate model fit diagnostics and term importance
library(mgcv)
library(Metrics)
library(stringr)
library(pracma)
library(MuMIn)

modFam = tw
masterDF = data.frame(read.csv("MasterWeeklyDF_plusJAX.csv"))
specs = cbind(c("SBCD","Risso","SFPW","Blainville","Gervais","Cuvier","Sowerby","True","Kogia","SpermWhale"),
              c("Dd","Gg","Gm","Md","Me","Zc","Mb","Mm","Kg","Pm"))

# Center and scale all predictors
masterDF[,3:12] = scale(masterDF[,3:12],center=TRUE,scale=TRUE)

# Round presence data
masterDF[,13:22] = round(masterDF[,13:22])

# Remove incomplete observations (NAs in FSLE)
badRows = unique(which(is.na(masterDF),arr.ind=TRUE)[,1])
masterDF = masterDF[-badRows,]

DE = matrix(ncol=1,nrow=dim(specs)[1])
colnames(DE) = "DE"
rownames(DE) = specs[,2]
AllModStats = matrix(ncol=6,nrow=dim(specs)[1])
colnames(AllModStats) = c("Rho","MAE","MAENorm","MAPE","RMSE","RMSENorm")
rownames(AllModStats) = specs[,2]
NonZeroModStats = matrix(ncol=6,nrow=dim(specs)[1])
colnames(NonZeroModStats) = c("Rho","MAE","MAENorm","MAPE","RMSE","RMSENorm")
rownames(NonZeroModStats) = specs[,2]

for (i in 1:dim(specs)[1]){
  
  # get species abbreviation
  abbrev = specs[i,2]
  
  # load model(s)
  load(paste(getwd(),'/',specs[i,1],'_WeeklyRegionalModel_Updated.Rdata',sep=""))
  
  # calculate deviance explained, averaging across models if necessary
  if (numel(names(topMods))>1){
    DE_list = numeric()
    AIC_list = numeric()
    for (j in 1:numel(names(topMods))){
      DE_list = rbind(DE_list,1-(topMods[[j]]$deviance/topMods[[j]]$null.deviance))
      AIC_list = rbind(AIC_list,topMods[[j]]$aic)
    }
    DE[i] = weighted.mean(DE_list,-AIC_list)
  }else{DE[i] = 1-(topMods[[1]]$deviance/topMods[[1]]$null.deviance)}
  
  # Retrain model with 2/3 of data
  trainInd = sample(1:nrow(masterDF),floor(nrow(masterDF)*.66))
  testInd = setdiff(1:nrow(masterDF),trainInd)
  masterDF$Pres = masterDF[[abbrev]]
  
  if (numel(names(topMods))>1){ # average top models, if necessary
    valModList = list()
    for (j in 1:numel(names(topMods))){
      valMod = update(topMods[[j]],data=masterDF[trainInd,])
      valModList[[j]] = valMod
    }
    valMod = model.avg(valModList,subset=delta<2,fit=TRUE)
  }else{valMod = update(optWeekMod,data=masterDF[trainInd,])}
  
  # predict on validation data
  modPreds = predict(valMod,masterDF[testInd,],full=TRUE,type="response",backtransform=FALSE)
  
  # get true values for validation data
  true = masterDF$Pres[testInd]

  # find non-zero validation data to compare to associated fitted values
  pres = which(masterDF$Pres[testInd]>0)
  
  # calculate Spearman's rank correlation
  NonZeroModStats[i,1] = round(cor(masterDF$Pres[testInd[pres]], modPreds[pres],method="spearman"),digits=4)
  AllModStats[i,1] = round(cor(masterDF$Pres[testInd], modPreds,method="spearman"),digits=4)
  
  # calculate mean absolute error
  NonZeroModStats[i,2] = round(mean(abs(masterDF$Pres[testInd[pres]]- modPreds[pres])),digits=4)
  AllModStats[i,2] = round(mean(abs(masterDF$Pres[testInd]- modPreds)),digits=4)
  
  # normalize MAE by 10th-90th percentile range
  quants = quantile(masterDF$Pres[testInd[pres]],probs=c(0.1,0.9))
  Ziqr = quants[2]-quants[1]
  NonZeroModStats[i,3] = round(NonZeroModStats[i,2]/Ziqr,digits=4)
  quants = quantile(masterDF$Pres[testInd],probs=c(0.1,0.9))
  Aiqr = quants[2]-quants[1]
  AllModStats[i,3] = round(AllModStats[i,2]/Aiqr,digits=4)

  # calculate mean absolute percent error
  NonZeroModStats[i,4] = round(mean(abs(masterDF$Pres[testInd[pres]]- modPreds[pres])/masterDF$Pres[testInd[pres]]),digits=4)
  AllModStats[i,4] = round(mean(abs(masterDF$Pres[testInd]- modPreds)/masterDF$Pres[testInd]),digits=4) # get -Inf from dividing by 0

  # calculate root mean square error
  NonZeroModStats[i,5] = round(mean((masterDF$Pres[testInd[pres]]-modPreds[pres])^2),digits=4)
  AllModStats[i,5] = round(mean((masterDF$Pres[testInd]-modPreds)^2),digits=4)

  # normalize RMSE by 10th-90th percentile range
  NonZeroModStats[i,6] = round(NonZeroModStats[i,5]/Ziqr,digits=4)
  AllModStats[i,6] = round(AllModStats[i,5]/Aiqr,digits=4)
  
  # # re-plot model
  #   thisForm = as.character(optWeekMod$formula)[3]
  #   startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
  #   termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
  #   termInd = c(0,termInd,str_length(thisForm)+1)
  #   allTerms = character()
  #   for (j in 1:length(termInd)-1){
  #     thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
  #     allTerms = c(allTerms,thisTerm)
  #   }
  #   nullTerm = str_which(allTerms,"1")
  #   allTerms = allTerms[-nullTerm]
  #   if (length(allTerms)<=4){
  #     wd = 600
  #   } else if (length(allTerms)>=5){
  #     wd=900
  #   }
  # png(filename=paste(specs[i,1],'_allSitesWeekly.png',sep=""),width=wd,height=600)
  # plot.gam(optWeekMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
  # while (dev.cur()>1) {dev.off()}
  # pdf(file=paste(specs[i,1],'_allSitesWeekly_forPrint.pdf',sep=""),width=wd/100*0.75,height=4.5,pointsize=10)
  # plot.gam(optWeekMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
  # while (dev.cur()>1) {dev.off()}
}

ModStats = cbind(NonZeroModStats,AllModStats,DE)

# Save fit statistics
write.csv(ModStats,'ModelStatistics_Updated.csv',row.names=TRUE)

#   # Calculate deviance explained by each model term -----------
# 
#   # model deviance
#   modDev = optWeekMod$deviance
# 
#   # null deviance
#   nullDev = optWeekMod$null.deviance
# 
#   # get model terms
#   thisForm = as.character(optWeekMod$formula)[3]
#   startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
#   termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
#   termInd = c(0,termInd,str_length(thisForm)+1)
#   allTerms = character()
#   for (j in 1:length(termInd)-1){
#     thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
#     allTerms = c(allTerms,thisTerm)
#   }
#   
#   devExp = matrix(nrow=length(allTerms),ncol=1)
# 
#   for (k in 1:length(allTerms)){
#     newForm = as.formula(paste(specs[i,2]," ~ ",paste(allTerms[-k],collapse="+"),sep=""))
#     
#     redMod<-gam(newForm,
#                 data=masterDF,
#                 family=tw,
#                 gamma=1.4,
#                 na.action="na.fail",
#                 method="REML",
#                 select=TRUE,
#                 sp=optWeekMod$sp[-k])
#     redDev = redMod$deviance
#     devExp2[k] = ((redDev-modDev)/(nullDev))
#   }
#   
#   termConts = data.frame(allTerms,devExp)
#   termConts$prop = termConts$devExp/sum(termConts$devExp)
#   termConts$cont = (termConts$prop*summary(optWeekMod)$dev.expl)*100
#   
#   rownames(devExp)=allTerms
#  devExp = cbind(devExp,NA,NA)
#  devExp[,2] = devExp[,1]/sum(devExp[,1])
#  devExp[,3] = devExp[,2]*summary(optWeekMod)$dev.expl*100  
#  
# }
