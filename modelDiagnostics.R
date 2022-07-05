# Calculate model fit diagnostics and term importance
library(mgcv)
library(Metrics)
library(stringr)
modFam = tw
masterDF = data.frame(read.csv("MasterWeeklyDF.csv"))
specs = cbind(rev(c("SBCD","Risso","SFPW","Blainville","Gervais","Cuvier","Sowerby","True","Kogia","SpermWhale")),
              rev(c("Dd","Gg","Gm","Md","Me","Zc","Mb","Mm","Kg","Pm")))


# Remove zeros in FSLE data to prepare for later transformation
masterDF$FSLE0[masterDF$FSLE0==0] = NA

# Transform data to fix skew, get all predictors on a similar scale
masterDF$log_Chl0 = log10(masterDF$Chl0)
masterDF$log_abs_FSLE0 = log10(abs(masterDF$FSLE0))
masterDF$sqrt_CEddyDist0 = sqrt(masterDF$CEddyDist0)
masterDF$sqrt_AEddyDist0 = sqrt(masterDF$AEddyDist0)
masterDF$sqrt_VelAsp0 = sqrt(masterDF$VelAsp0)
# masterDF$sqrt_VelAsp700 = sqrt(masterDF$VelAsp700)
masterDF$sqrt_EKE0 = sqrt(masterDF$EKE0)

# Remove incomplete observations (NAs in FSLE)
badRows = unique(which(is.na(masterDF),arr.ind=TRUE)[,1])
masterDF = masterDF[-badRows,]

ModStats = matrix(nrow=4,ncol=dim(specs)[1])
rownames(ModStats) = c("Rho","MAE","Norm_MAE","MAPE")
colnames(ModStats) = specs[,2]

for (i in 1:dim(specs)[1]){
  # load model
  tempMod = load(paste(getwd(),'/',specs[i,1],'_WeeklyRegionalModel.Rdata',sep=""))
  
  # find non-zero input data to compare to associated fitted values
  abbrev = specs[i,2]
  pres = which(masterDF[[abbrev]]>0)
  
  # calculate model predictions for input data
  modPreds = predict.gam(optWeekMod,masterDF,type="response") # should this be method="response"??
  
  # calculate Spearman's rank correlation
  ModStats[1,i] = cor(masterDF[[abbrev]][pres], modPreds[pres],method="spearman")
  
  # calculate mean absolute error
  ModStats[2,i] = mean(abs(masterDF[[abbrev]][pres]- modPreds[pres]))
  
  # normalize MAE
  ModStats[3,i] = ModStats[2,i]/max(masterDF[[abbrev]][pres])

  # calculate mean absolute percent error
  ModStats[4,i] = mean(abs((masterDF[[abbrev]][pres]-modPreds[pres])/masterDF[[abbrev]][pres]))
}  
  # Calculate deviance explained by each model term -----------
  
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
#   devExp = matrix(nrow=length(allTerms),ncol=1)
#   
#   for (i in 1:length(allTerms)){
#     redMod<-eval(parse(text=paste("update(optWeekMod, . ~ . - ", allTerms[i], ")", sep="")))
#     redDev = redMod$deviance
#     devExp[i] = ((redDev-modDev)/(nullDev))*100
#   }
#   
#   rownames(devExp)=c("FSLE0","VelAsp0","Temp0","Temp700","Sal0","Intercept")
# }
