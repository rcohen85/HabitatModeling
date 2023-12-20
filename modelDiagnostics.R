# Calculate model fit diagnostics and term importance
library(mgcv)
library(Metrics)
library(stringr)
modFam = tw
masterDF = data.frame(read.csv("MasterWeeklyDF.csv"))
specs = cbind(c("SBCD","Risso","SFPW","Blainville","Gervais","Cuvier","Sowerby","True","Kogia","SpermWhale"),
              c("Dd","Gg","Gm","Md","Me","Zc","Mb","Mm","Kg","Pm"))


# # Remove zeros in FSLE data to prepare for later transformation
# masterDF$FSLE0[masterDF$FSLE0==0] = NA
# 
# # Transform data to fix skew, get all predictors on a similar scale
# masterDF$log_Chl0 = log10(masterDF$Chl0)
# masterDF$log_abs_FSLE0 = log10(abs(masterDF$FSLE0))
# masterDF$sqrt_CEddyDist0 = sqrt(masterDF$CEddyDist0)
# masterDF$sqrt_AEddyDist0 = sqrt(masterDF$AEddyDist0)
# masterDF$sqrt_VelAsp0 = sqrt(masterDF$VelAsp0)
# masterDF$sqrt_EKE0 = sqrt(masterDF$EKE0)

# Center and scale all predictors
masterDF[,3:12] = scale(masterDF[,3:12],center=TRUE,scale=TRUE)

# Remove incomplete observations (NAs in FSLE)
badRows = unique(which(is.na(masterDF),arr.ind=TRUE)[,1])
masterDF = masterDF[-badRows,]
weeklyDF = masterDF

AllModStats = matrix(ncol=6,nrow=dim(specs)[1])
colnames(AllModStats) = c("Rho","MAE","MAENorm","MAPE","RMSE","RMSENorm")
rownames(AllModStats) = specs[,2]
ZeroModStats = matrix(ncol=6,nrow=dim(specs)[1])
colnames(ZeroModStats) = c("Rho","MAE","MAENorm","MAPE","RMSE","RMSENorm")
rownames(ZeroModStats) = specs[,2]

for (i in 1:dim(specs)[1]){
  # load model
  load(paste(getwd(),'/',specs[i,1],'_WeeklyRegionalModel.Rdata',sep=""))
  
  # find non-zero input data to compare to associated fitted values
  abbrev = specs[i,2]
  pres = which(masterDF[[abbrev]]>0)
  
  # calculate model predictions for input data
  modPreds = predict.gam(optWeekMod,masterDF,type="response")
  
  # calculate Spearman's rank correlation
  ZeroModStats[i,1] = round(cor(masterDF[[abbrev]][pres], modPreds[pres],method="spearman"),digits=4)
  AllModStats[i,1] = round(cor(masterDF[[abbrev]], modPreds,method="spearman"),digits=4)
  
  # calculate mean absolute error
  ZeroModStats[i,2] = round(mean(abs(masterDF[[abbrev]][pres]- modPreds[pres])),digits=4)
  AllModStats[i,2] = round(mean(abs(masterDF[[abbrev]]- modPreds)),digits=4)
  
  # normalize MAE by 10th-90th percentile range
  quants = quantile(masterDF[[abbrev]][pres],probs=c(0.1,0.9))
  Ziqr = quants[2]-quants[1]
  ZeroModStats[i,3] = round(ZeroModStats[i,2]/Ziqr,digits=4)
  quants = quantile(masterDF[[abbrev]],probs=c(0.1,0.9))
  Aiqr = quants[2]-quants[1]
  AllModStats[i,3] = round(AllModStats[i,2]/Aiqr,digits=4)

  # calculate mean absolute percent error
  ZeroModStats[i,4] = round(mean(abs(masterDF[[abbrev]][pres]- modPreds[pres])/masterDF[[abbrev]][pres]),digits=4)
  AllModStats[i,4] = round(mean(abs(masterDF[[abbrev]]- modPreds)/masterDF[[abbrev]]),digits=4) # get -Inf from dividing by 0

  # calculate root mean square error
  ZeroModStats[i,5] = round(mean((masterDF[[abbrev]][pres]-modPreds[pres])^2),digits=4)
  AllModStats[i,5] = round(mean((masterDF[[abbrev]]-modPreds)^2),digits=4)

  # calculate root mean square error
  ZeroModStats[i,6] = round(ZeroModStats[i,5]/Ziqr,digits=4)
  AllModStats[i,6] = round(AllModStats[i,5]/Aiqr,digits=4)
  
  # re-plot model
    thisForm = as.character(optWeekMod$formula)[3]
    startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
    termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
    termInd = c(0,termInd,str_length(thisForm)+1)
    allTerms = character()
    for (j in 1:length(termInd)-1){
      thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
      allTerms = c(allTerms,thisTerm)
    }
    nullTerm = str_which(allTerms,"1")
    allTerms = allTerms[-nullTerm]
    if (length(allTerms)<=4){
      wd = 600
    } else if (length(allTerms)>=5){
      wd=900
    }
  png(filename=paste(specs[i,1],'_allSitesWeekly.png',sep=""),width=wd,height=600)
  plot.gam(optWeekMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
  while (dev.cur()>1) {dev.off()}
  pdf(file=paste(specs[i,1],'_allSitesWeekly_forPrint.pdf',sep=""),width=wd/100*0.75,height=4.5,pointsize=10)
  plot.gam(optWeekMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
  while (dev.cur()>1) {dev.off()}
}

ModStats = cbind(ZeroModStats,AllModStats)

# Save fit statistics
write.csv(ModStats,'ModelStatistics.csv',row.names=TRUE)

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
