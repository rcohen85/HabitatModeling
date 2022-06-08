library(tidyverse)
# library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

## GAM approach ---------------------
# Regional model
spec = 'Kogia'
outDir = "E:/ModelingCovarData/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/Kogia_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

# create weekly time series to reduce autocorrelation
stDt = as.Date("2016-05-01")
edDt = as.Date("2019-04-30")
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS')
allDates = stDt:edDt
weekDates = seq.Date(stDt,edDt,by=7)
weeklyDF = as.numeric()

for (l in 1:length(sites)) {
  
  # Create dataframe to hold data (or NAs) for all dates
  fullDatesDF = data.frame(matrix(nrow=length(allDates), ncol=44))
  fullDatesDF[,1] = allDates
  
  # sort the observations we have for this site into the full date range
  thisSite = which(!is.na(str_match(data$Site,sites[l])))
  matchRow = match(data$Date[thisSite],allDates)
  fullDatesDF[matchRow,2:dim(data)[2]] = data[thisSite,2:dim(data)[2]]
  
  colnames(fullDatesDF) = colnames(data)
  
  # create grouping variable
  weekID = rep(1:length(weekDates),each=7)
  weekID = weekID[1:dim(fullDatesDF)[1]]
  fullDatesDF$WeekID = weekID
  
  summaryData = fullDatesDF %>%
    group_by(WeekID) %>%
    summarize(Pres=sum(Pres,na.rm=TRUE))
  
  summaryData$Site = sites[l]
  summaryData$WeekDate = weekDates
  
  for (j in 4:(dim(data)[2])){
    var = names(data)[j]
    # calculate weekly average for this covar
    eval(parse(text=paste('thisCovar=fullDatesDF%>%group_by(WeekID)%>%summarize(',var,'=mean(',var,',na.rm=TRUE))',sep="")))
    eval(parse(text='summaryData[[var]]=unlist(thisCovar[,2])'))
    
  }
  weeklyDF = rbind(weeklyDF,summaryData)
}

# Remove zeros in FSLE data to prepare for later transformation
data$FSLE0[data$FSLE0==0] = NA
weeklyDF$FSLE0[weeklyDF$FSLE0==0] = NA

# Remove incomplete observations (NAs in FSLE)
badRows = which(is.na(data),arr.ind=TRUE)[,1]
data = data[-badRows,]
badRows = unique(which(is.na(weeklyDF),arr.ind=TRUE)[,1])
weeklyDF = weeklyDF[-badRows,]


# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c ("Chl0",
                   "EKE0",
                   "FSLE0",
                   "Sal0",
                   "Sal700",
                   "SSH0",
                   "Temp0",
                   "Temp700",
                   "VelAsp0",
                   "VelMag0",
                   "CEddyDist0")


# determine whether covars should be included as linear or smooth terms
# seafloor aspect and water velocity direction are included as cyclic smooths

modOpts = c("linMod","threeKnots","fourKnots","fiveKnots")
AIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"Asp")){
    bs = "cc"
  } else { bs = "cs"}
  
  modelCall = paste('gam(Pres~data$',smoothVarList[i],',data=data,family=poisson)',sep="")
  linMod = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=3),data=data,family=poisson)',sep="")
  smoothMod1 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=4),data=data,family=poisson)',sep="")
  smoothMod2 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=5),data=data,family=poisson)',sep="")
  smoothMod3 = eval(parse(text=modelCall))
  
  AIC_votes[i,1:4] = c(AIC(linMod)[[1]],AIC(smoothMod1)[[1]],AIC(smoothMod2)[[1]],AIC(smoothMod3)[[1]])
  AIC_votes[i,5] = modOpts[which.min(AIC_votes[i,1:4])]
}

colnames(AIC_votes) = c(modOpts,"Best")
rownames(AIC_votes) = smoothVarList[]
AIC_votes

#               linMod             threeKnots         fourKnots          fiveKnots          Best       
# Chl0       "14357.1508790267" "14187.8242065542" "14187.2171494088" "13941.1524678567" "fiveKnots"
# EKE0       "16638.7449683075" "16521.4655035031" "16378.8487351627" "16368.6123892055" "fiveKnots"
# FSLE0      "16875.4204340965" "16699.8589647431" "16696.3977011661" "16694.2813553344" "fiveKnots"
# Sal0       "14035.1122567426" "13711.8024362691" "13662.2162646885" "13652.3242698491" "fiveKnots"
# Sal700     "13900.6377156938" "13703.2405404535" "13683.0183061573" "13668.9222207876" "fiveKnots"
# SSH0       "13724.0035388473" "13490.5156793078" "13333.7463268739" "13319.2132463421" "fiveKnots"
# Temp0      "15339.8502379472" "15161.3724603807" "15052.8197607025" "15069.2150725657" "fourKnots"
# Temp700    "14436.5705360095" "14216.9153425856" "13939.8221224786" "13856.516560931"  "fiveKnots"
# VelAsp0    "16653.4613165684" "16387.3045168356" "16387.3045168356" "16372.7059769484" "fiveKnots"
# VelMag0    "16935.9039497736" "16700.667914537"  "16641.2919299618" "16624.8167730333" "fiveKnots"
# CEddyDist0 "15917.3301485828" "15916.540055656"  "15893.3462341004" "15892.5226047707" "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
              + s(log(Chl0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4)
              + s(Temp700,bs="cs",k=5)
              + s(sqrt(VelAsp0),bs="cs",k=5)
              + s(VelMag0,bs="cs",k=5),
              data=data,
              family=poisson,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
              + s(log(Chl0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              + s(Sal0,bs="cs",k=5)
              # + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4)
              + s(Temp700,bs="cs",k=5)
              + s(sqrt(VelAsp0),bs="cs",k=5)
              + s(VelMag0,bs="cs",k=5),
              data=weeklyDF,
              family=poisson,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

# check convergence
fullMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                     para s(sqrt(CEddyDist0)) s(log(Chl0)) s(log(abs(FSLE0))) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(sqrt(VelAsp0)) s(VelMag0)
# para                   1              0.0000       0.0000             0.0000  0.0000    0.0000  0.0000   0.0000     0.0000           0.0000     0.0000
# s(sqrt(CEddyDist0))    0              1.0000       0.0710             0.0339  0.1544    0.1578  0.1012   0.1079     0.1444           0.0652     0.1003
# s(log(Chl0))           0              0.0772       1.0000             0.0429  0.2279    0.2208  0.2121   0.4680     0.2419           0.0798     0.0391
# s(log(abs(FSLE0)))     0              0.0403       0.0481             1.0000  0.1152    0.1119  0.1212   0.0858     0.1250           0.0878     0.1094
# s(Sal0)                0              0.1665       0.2707             0.1061  1.0000    0.6566  0.2936   0.4170     0.5783           0.2216     0.2566
# s(Sal700)              0              0.1669       0.2947             0.1096  0.6580    1.0000  0.3471   0.3848     0.6630           0.2269     0.2893
# s(SSH0)                0              0.1550       0.3703             0.1005  0.3984    0.3781  1.0000   0.4919     0.4560           0.2213     0.3099
# s(Temp0)               0              0.0961       0.3702             0.0549  0.3188    0.3165  0.2625   1.0000     0.3924           0.1228     0.1350
# s(Temp700)             0              0.1620       0.2989             0.1076  0.5665    0.6750  0.3340   0.4326     1.0000           0.2122     0.1966
# s(sqrt(VelAsp0))       0              0.0628       0.0754             0.0722  0.1989    0.1884  0.1666   0.1491     0.1923           1.0000     0.2646
# s(VelMag0)             0              0.0769       0.0322             0.0522  0.1606    0.1537  0.2403   0.1091     0.1375           0.1721     1.0000

# Sal0 concurved with Sal700 and Temp700
# Sal700 concurved with Sal0 and Temp700
# Temp700 concurved with Sal0 and Sal700
# check pvals of Sal0, Sal700, Temp700

# check p vals of concurved covars
format.pval(summary(gam(Pres~s(Sal0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Sal700,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Temp700,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16

# all p vlaues are the same :(
# remove sal0, temp700


# run reduced model
dayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
             + s(log(Chl0),bs="cs",k=5)
             + s(log(abs(FSLE0)),bs="cs",k=5)
             # + s(Sal0,bs="cs",k=5)
             + s(Sal700,bs="cs",k=5)
             + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=4)
             # + s(Temp700,bs="cs",k=5)
             + s(sqrt(VelAsp0),bs="cs",k=5)
             + s(VelMag0,bs="cs",k=5),
             data=data,
             family=poisson,
             method="REML",
             select=TRUE,
             gamma=1.4,
             na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
              + s(log(Chl0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              # + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4)
              # + s(Temp700,bs="cs",k=5)
              + s(sqrt(VelAsp0),bs="cs",k=5)
              + s(VelMag0,bs="cs",k=5),
              data=weeklyDF,
              family=poisson,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

# check convergence
dayMod$converged
# TRUE
weekMod$converged
# TRUE

conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                       para s(sqrt(CEddyDist0)) s(log(Chl0)) s(log(abs(FSLE0))) s(Sal700) s(SSH0) s(Temp0) s(sqrt(VelAsp0)) s(VelMag0)
# para                   1              0.0000       0.0000             0.0000    0.0000  0.0000   0.0000           0.0000     0.0000
# s(sqrt(CEddyDist0))    0              1.0000       0.0843             0.0554    0.1902  0.1092   0.1247           0.0701     0.1179
# s(log(Chl0))           0              0.0942       1.0000             0.0723    0.2177  0.2120   0.4581           0.0835     0.0758
# s(log(abs(FSLE0)))     0              0.0671       0.0759             1.0000    0.1727  0.1837   0.1155           0.1099     0.1353
# s(Sal700)              0              0.2016       0.3001             0.1782    1.0000  0.4132   0.3853           0.2162     0.3869
# s(SSH0)                0              0.1840       0.3717             0.1607    0.4023  1.0000   0.4799           0.2060     0.3791
# s(Temp0)               0              0.1156       0.3568             0.0859    0.3104  0.2578   1.0000           0.1111     0.1776
# s(sqrt(VelAsp0))       0              0.1123       0.1054             0.1417    0.2736  0.2668   0.1896           1.0000     0.4542
# s(VelMag0)             0              0.1072       0.0637             0.0882    0.2204  0.3237   0.1374           0.2034     1.0000

# everything looks good! all values below 0.5

dayModCompTable = dredge(dayMod,
                         beta="none",
                         evaluate=TRUE,
                         trace=TRUE)

weekModCompTable = dredge(weekMod,
                          beta="none",
                          evaluate=TRUE,
                          trace=TRUE)

# run optimal models
optDayMod = get.models(dayModCompTable,subset=1)
optDayMod = optDayMod[[names(optDayMod)]]
save(optDayMod,dayModCompTable,file=paste(outDir,'/',spec,'/','DailyRegionalModel.Rdata',sep=""))
optWeekMod = get.models(weekModCompTable,subset=1)
optWeekMod = optWeekMod[[names(optWeekMod)]]
save(optWeekMod,weekModCompTable,file=paste(outDir,'/',spec,'/','WeeklyRegionalModel.Rdata',sep=""))

# check p-values
summary(optDayMod)

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(log(Chl0), 
#                                                   bs = "cs", k = 5) + s(Sal700, bs = "cs", k = 5) + 
#   s(sqrt(CEddyDist0), bs = "cs", k = 5) + s(sqrt(VelAsp0), 
#                                             bs = "cs", k = 5) + s(SSH0, bs = "cs", k = 5) + 
#   s(Temp0, bs = "cs", k = 4) + s(VelMag0, bs = "cs", 
#                                  k = 5) + 1
# 
# Parametric coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.12169    0.03812  -55.66   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                        edf Ref.df  Chi.sq  p-value    
#   s(log(abs(FSLE0)))  3.1820      4  27.931 1.33e-06 ***
#   s(log(Chl0))        3.8406      4 123.049  < 2e-16 ***
#   s(Sal700)           0.9284      4   4.282   0.0149 *  
#   s(sqrt(CEddyDist0)) 2.9235      4  55.528  < 2e-16 ***
#   s(sqrt(VelAsp0))    2.2423      4   8.063   0.0146 *  
#   s(SSH0)             3.2541      4 221.357  < 2e-16 ***
#   s(Temp0)            2.7327      3 156.350  < 2e-16 ***
#   s(VelMag0)          3.1196      4  44.798  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.0794   Deviance explained = 31.1%
# -REML = 4623.8  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(log(Chl0), 
#                                                   bs = "cs", k = 5) + s(sqrt(CEddyDist0), bs = "cs", 
#                                                                         k = 5) + s(sqrt(VelAsp0), bs = "cs", k = 5) + s(SSH0, 
#                                                                                                                         bs = "cs", k = 5) + s(Temp0, bs = "cs", k = 4) + 
#   s(VelMag0, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.20792    0.03821  -5.442 5.26e-08 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                         edf Ref.df  Chi.sq  p-value    
#   s(log(abs(FSLE0)))  2.187      4   6.953   0.0226 *  
#   s(log(Chl0))        3.744      4  98.424  < 2e-16 ***
#   s(sqrt(CEddyDist0)) 2.963      4  87.596  < 2e-16 ***
#   s(sqrt(VelAsp0))    2.048      4  24.199 2.15e-06 ***
#   s(SSH0)             3.504      4 457.389  < 2e-16 ***
#   s(Temp0)            2.799      3 169.372  < 2e-16 ***
#   s(VelMag0)          3.364      4 108.882  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.299   Deviance explained =   54%
# -REML = 2042.5  Scale est. = 1         n = 1509

# plot
png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesDaily.png',sep=""),width=600,height=600)
plot.gam(optDayMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
while (dev.cur()>1) {dev.off()}
png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesWeekly.png',sep=""),width=600,height=600)
plot.gam(optWeekMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
while (dev.cur()>1) {dev.off()}


# Site-specific models ---------------
sites = unique(data$Site)
siteDayModList = list()
pValDayList = list()
siteDayModCompList = list()
siteWeekModList = list()
pValWeekList = list()
siteWeekModCompList = list()
for (i in 1:length(sites)){
  
  dayInd = which(!is.na(str_match(data$Site,sites[i])))
  dayData = data[dayInd,]
  weekInd = which(!is.na(str_match(weeklyDF$Site,sites[i])))
  weekData = weeklyDF[weekInd,]
  
  if (sum(dayData$Pres>0)>25){
    
    # run full daily model for this site
    fullSiteDayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
                         + s(log(Chl0),bs="cs",k=5)
                         + s(log(abs(FSLE0)),bs="cs",k=5)
                         # + s(Sal0,bs="cs",k=5)
                         + s(Sal700,bs="cs",k=5)
                         + s(SSH0,bs="cs",k=5)
                         + s(Temp0,bs="cs",k=4)
                         # + s(Temp700,bs="cs",k=5)
                         + s(sqrt(VelAsp0),bs="cs",k=5)
                         + s(VelMag0,bs="cs",k=5),
                         data=dayData,
                         family=poisson,
                         method="REML",
                         select=TRUE,
                         gamma=1.4,
                         na.action="na.fail")
    
    siteDayModCompTable = dredge(fullSiteDayMod,
                                 beta="none",
                                 evaluate=TRUE,
                                 trace=TRUE)
    siteDayModCompList[[sites[i]]] = siteDayModCompTable
    
    optSiteDayMod = get.models(siteDayModCompTable,subset=1)
    optSiteDayMod = optSiteDayMod[[names(optSiteDayMod)]]
    siteDayPV = summary(optSiteDayMod)$s.pv
    
    if (any(siteDayPV>=0.05)){ # Remove non-significant terms & re-run model iteratively until only signif covars remain
      flag = 1
      while (flag==1){
        # get terms from formula as strings
        thisForm = as.character(optSiteDayMod$formula)[3]
        startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
        termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
        termInd = c(0,termInd,str_length(thisForm)+1)
        allTerms = character()
        for (j in 1:length(termInd)-1){
          thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
          allTerms = c(allTerms,thisTerm)
        }
        # identify which terms were non-significant
        badVars = allTerms[siteDayPV>=0.05]
        dontNeed = which(!is.na(str_match(badVars,"1")))
        if (!is_empty(dontNeed)){
          badVars = badVars[-dontNeed]}
        # update model
        optSiteDayMod<-eval(parse(text=paste("update(optSiteDayMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
        siteDayPV = summary(optSiteDayMod)$s.pv
        if (!any(siteDayPV>=0.05)){
          siteDayModList[[sites[i]]] = optSiteDayMod
          pValDayList[[sites[i]]] = siteDayPV
          flag=0
        }
      }
    } else {
      siteDayModList[[sites[i]]] = optSiteDayMod
      pValDayList[[sites[i]]] = siteDayPV
    }
    
    sink(paste(outDir,'/',spec,'/',spec,'_',sites[i],'_DailySummary.txt',sep=""))
    print(summary(siteDayModList[[sites[i]]]))
    sink()
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'Daily.png',sep=""),width=600,height=600)
    plot.gam(siteDayModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
    while (dev.cur()>1) {dev.off()}
    
  }
  
  if (sum(weekData$Pres>0)>25){
    #run full weekly model for this site
    fullSiteWeekMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
                          + s(log(Chl0),bs="cs",k=5)
                          + s(log(abs(FSLE0)),bs="cs",k=5)
                          # + s(Sal0,bs="cs",k=5)
                          + s(Sal700,bs="cs",k=5)
                          + s(SSH0,bs="cs",k=5)
                          + s(Temp0,bs="cs",k=4)
                          # + s(Temp700,bs="cs",k=5)
                          + s(sqrt(VelAsp0),bs="cs",k=5)
                          + s(VelMag0,bs="cs",k=5),
                          data=weekData,
                          family=poisson,
                          method="REML",
                          select=TRUE,
                          gamma=1.4,
                          na.action="na.fail")
    
    siteWeekModCompTable = dredge(fullSiteWeekMod,
                                  beta="none",
                                  evaluate=TRUE,
                                  trace=TRUE)
    siteWeekModCompList[[sites[i]]] = siteWeekModCompTable
    
    optSiteWeekMod = get.models(siteWeekModCompTable,subset=1)
    optSiteWeekMod = optSiteWeekMod[[names(optSiteWeekMod)]]
    siteWeekPV = summary(optSiteWeekMod)$s.pv
    
    if (any(siteWeekPV>=0.05)){ # Remove non-significant terms & re-run model iteratively until only signif covars remain
      flag = 1
      while (flag==1){
        # get terms from formula as strings
        thisForm = as.character(optSiteWeekMod$formula)[3]
        startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
        termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
        termInd = c(0,termInd,str_length(thisForm)+1)
        allTerms = character()
        for (j in 1:length(termInd)-1){
          thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
          allTerms = c(allTerms,thisTerm)
        }
        # identify which terms were non-significant
        badVars = allTerms[siteWeekPV>=0.05]
        dontNeed = which(!is.na(str_match(badVars,"1")))
        if (!is_empty(dontNeed)){
          badVars = badVars[-dontNeed]}
        # update model
        optSiteWeekMod<-eval(parse(text=paste("update(optSiteWeekMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
        siteWeekPV = summary(optSiteWeekMod)$s.pv
        if (!any(siteWeekPV>=0.05)){
          siteWeekModList[[sites[i]]] = optSiteWeekMod
          pValWeekList[[sites[i]]] = siteWeekPV
          flag=0
        }
      }
    } else {
      siteWeekModList[[sites[i]]] = optSiteWeekMod
      pValWeekList[[sites[i]]] = siteWeekPV
    }
    
    sink(paste(outDir,'/',spec,'/',spec,'_',sites[i],'_WeeklySummary.txt',sep=""))
    print(summary(siteWeekModList[[sites[i]]]))
    sink()
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'Weekly.png',sep=""),width=600,height=600)
    plot.gam(siteWeekModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
    while (dev.cur()>1) {dev.off()}
  }
}

save(siteDayModList,pValDayList,siteDayModCompList,file=paste(outDir,'/',spec,'/','DailySiteSpecificModels.Rdata',sep=""))
save(siteWeekModList,pValWeekList,siteWeekModCompList,file=paste(outDir,'/',spec,'/','WeeklySiteSpecificModels.Rdata',sep=""))
