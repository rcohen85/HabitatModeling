library(tidyverse)
library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)
library(forecast)
library(nlme)
library(itsadug)

## GAM approach ---------------------
# Regional model
spec = 'Blainville'
outDir = "J:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('J:/Chpt_3/ModelData/Blainville_masterDF.csv'))
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
smoothVarList = c ("EKE0",
                   "Sal0",
                   "Sal700",
                   "Temp0",
                   "VelAsp700",
                   "VelMag700",
                   "AEddyDist0",
                   "CEddyDist0")

# Test for how a term should be included in the model
modOpts = c("linMod","threeKnots","fourKnots","fiveKnots")
AIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"Asp")){
    bs = "cc"
  } else { bs = "cs"}
  
  modelCall = paste('gam(Pres~data$',smoothVarList[i],',data=data,family=poisson,method="REML",select=TRUE)',sep="")
  linMod = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=3),data=data,family=poisson,method="REML",select=TRUE)',sep="")
  smoothMod1 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=4),data=data,family=poisson,method="REML",select=TRUE)',sep="")
  smoothMod2 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=5),data=data,family=poisson,method="REML",select=TRUE)',sep="")
  smoothMod3 = eval(parse(text=modelCall))
  
  AIC_votes[i,1:4] = c(AIC(linMod)[[1]],AIC(smoothMod1)[[1]],AIC(smoothMod2)[[1]],AIC(smoothMod3)[[1]])
  AIC_votes[i,5] = modOpts[which.min(AIC_votes[i,1:4])]
}

colnames(AIC_votes) = c(modOpts,"Best")
rownames(AIC_votes) = smoothVarList[]
AIC_votes

#               linMod             threeKnots         fourKnots          fiveKnots          Best        
# EKE0       "50229.7821840425" "49102.9411588034" "48841.8266393599" "47684.7291157758" "fiveKnots" 
# Sal0       "26741.5411590067" "26724.6003965425" "26725.7023282191" "26728.2247649327" "threeKnots"
# Sal700     "24825.5752374581" "21735.0388546467" "21735.0004685225" "21734.8412500601" "fiveKnots" 
# Temp0      "44737.1144208223" "41509.3188145839" "40525.3030201052" "40253.0119236911" "fiveKnots" 
# VelAsp700  "49951.691572375"  "49406.3876230599" "49406.3876230599" "45115.6288995681" "fiveKnots" 
# VelMag700  "48786.5628418564" "47805.8093871974" "47705.0533853665" "47679.6627447923" "fiveKnots" 
# AEddyDist0 "43582.4919291093" "41693.166901377"  "41213.6958439186" "41208.5188702519" "fiveKnots" 
# CEddyDist0 "47912.9011422381" "47864.315207453"  "47799.2114943758" "47793.8887525785" "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(sqrt(EKE0),bs="cs",k=5)
              + s(Sal0,bs="cs",k=3)
              + s(Sal700,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(sqrt(VelAsp700),bs="cc",k=5)
              + s(VelMag700,bs="cs",k=5)
              + s(sqrt(AEddyDist0),bs="cs",k=5)
              + s(sqrt(CEddyDist0),bs="cs",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              select=TRUE)

weekMod = gam(Pres ~ s(sqrt(EKE0),bs="cs",k=5)
              + s(Sal0,bs="cs",k=3)
              + s(Sal700,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(sqrt(VelAsp700),bs="cc",k=5)
              + s(VelMag700,bs="cs",k=5)
              + s(sqrt(AEddyDist0),bs="cs",k=5)
              + s(sqrt(CEddyDist0),bs="cs",k=5),
              data=weeklyDF,
              family=poisson,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              select=TRUE)

# check convergence
fullMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                       para s(sqrt(EKE0)) s(Sal0) s(Sal700) s(Temp0) s(sqrt(VelAsp700)) s(VelMag700) s(sqrt(AEddyDist0)) s(sqrt(CEddyDist0))
# para                   1        0.0000  0.0000    0.0000   0.0000             0.0000       0.0000              0.0000              0.0000
# s(sqrt(EKE0))          0        1.0000  0.3263    0.1889   0.0897             0.4977       0.3251              0.0387              0.0641
# s(Sal0)                0        0.1887  1.0000    0.5712   0.2023             0.2004       0.1948              0.0555              0.1567
# s(Sal700)              0        0.2259  0.8873    1.0000   0.2254             0.2350       0.3346              0.1036              0.1669
# s(Temp0)               0        0.1241  0.4350    0.3239   1.0000             0.1193       0.1582              0.0462              0.1003
# s(sqrt(VelAsp700))     0        0.5048  0.3814    0.2160   0.0961             1.0000       0.3186              0.0481              0.0766
# s(VelMag700)           0        0.2084  0.2724    0.1909   0.0882             0.1834       1.0000              0.0133              0.0905
# s(sqrt(AEddyDist0))    0        0.0402  0.0980    0.0598   0.0251             0.0480       0.0091              1.0000              0.0185
# s(sqrt(CEddyDist0))    0        0.0659  0.2737    0.1578   0.0673             0.0753       0.1171              0.0181              1.0000

# EKE0 problematic with VelAsp700
# Sal0 problematic w Sal700, Temp0, VelAsp700, EKE
# Sal700 problematic w Sal0, less w Temp0

format.pval(summary(gam(Pres~s(EKE0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# "< 2.220446e-16"
format.pval(summary(gam(Pres~s(Sal0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# "< 2.220446e-16"
format.pval(summary(gam(Pres~s(Sal700,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# "< 2.220446e-16"

# All equally significant, arbitrarily taking out Sal0, EKE0

# run reduced model
dayMod = gam(Pres ~ 
               # s(sqrt(EKE0),bs="cs",k=5)
             # + s(Sal0,bs="cs",k=3)
             + s(Sal700,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5)
             + s(sqrt(VelAsp700),bs="cc",k=5)
             + s(VelMag700,bs="cs",k=5)
             + s(sqrt(AEddyDist0),bs="cs",k=5)
             + s(sqrt(CEddyDist0),bs="cs",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              select=TRUE)

weekMod = gam(Pres ~ 
                # s(sqrt(EKE0),bs="cs",k=5)
                # + s(Sal0,bs="cs",k=3)
                + s(Sal700,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(sqrt(VelAsp700),bs="cc",k=5)
              + s(VelMag700,bs="cs",k=5)
              + s(sqrt(AEddyDist0),bs="cs",k=5)
              + s(sqrt(CEddyDist0),bs="cs",k=5),
             data=weeklyDF,
             family=poisson,
             gamma=1.4,
             na.action="na.fail",
             method="REML",
             select=TRUE)

# check convergence
dayMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(dayMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                       para s(Sal700) s(Temp0) s(sqrt(VelAsp700)) s(VelMag700) s(sqrt(AEddyDist0)) s(sqrt(CEddyDist0))
# para                   1    0.0000   0.0000             0.0000       0.0000              0.0000              0.0000
# s(Sal700)              0    1.0000   0.2254             0.2350       0.3346              0.1036              0.1669
# s(Temp0)               0    0.3239   1.0000             0.1193       0.1582              0.0462              0.1003
# s(sqrt(VelAsp700))     0    0.2160   0.0961             1.0000       0.3186              0.0481              0.0766
# s(VelMag700)           0    0.1909   0.0882             0.1834       1.0000              0.0133              0.0905
# s(sqrt(AEddyDist0))    0    0.0598   0.0251             0.0480       0.0091              1.0000              0.0185
# s(sqrt(CEddyDist0))    0    0.1578   0.0673             0.0753       0.1171              0.0181              1.0000

# only slight concurvity remains

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
#   Pres ~ s(Sal700, bs = "cs", k = 5) + s(sqrt(AEddyDist0), bs = "cs", 
#                                          k = 5) + s(sqrt(CEddyDist0), bs = "cs", k = 5) + s(sqrt(VelAsp700), 
#                                                                                             bs = "cc", k = 5) + s(Temp0, bs = "cs", k = 5) + s(VelMag700, 
#                                                                                                                                                bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -31.49      14.21  -2.216   0.0267 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(Sal700)           1.999      4 1267.4  <2e-16 ***
#   s(sqrt(AEddyDist0)) 3.895      4  511.8  <2e-16 ***
#   s(sqrt(CEddyDist0)) 2.657      4  175.4  <2e-16 ***
#   s(sqrt(VelAsp700))  2.983      3  960.6  <2e-16 ***
#   s(Temp0)            2.989      4  331.4  <2e-16 ***
#   s(VelMag700)        3.796      4  574.3  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.421   Deviance explained = 72.5%
# -REML = 5950.5  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(Sal700, bs = "cs", k = 5) + s(sqrt(AEddyDist0), bs = "cs", 
#                                          k = 5) + s(sqrt(CEddyDist0), bs = "cs", k = 5) + s(sqrt(VelAsp700), 
#                                                                                             bs = "cc", k = 5) + s(Temp0, bs = "cs", k = 5) + s(VelMag700, 
#                                                                                                                                                bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -52.57      25.54  -2.059   0.0395 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(Sal700)           1.996      4 1211.5  <2e-16 ***
#   s(sqrt(AEddyDist0)) 2.848      4  458.8  <2e-16 ***
#   s(sqrt(CEddyDist0)) 3.606      4  236.0  <2e-16 ***
#   s(sqrt(VelAsp700))  2.964      3  675.1  <2e-16 ***
#   s(Temp0)            2.950      4  261.1  <2e-16 ***
#   s(VelMag700)        3.845      4  560.0  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.608   Deviance explained = 82.2%
# -REML = 3129.3  Scale est. = 1         n = 1509

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
    fullSiteDayMod = gam(Pres ~ s(Sal700,bs="cs",k=5)
                         + s(Temp0,bs="cs",k=5)
                         + s(sqrt(VelAsp700),bs="cc",k=5)
                         + s(VelMag700,bs="cs",k=5)
                         + s(sqrt(AEddyDist0),bs="cs",k=5)
                         + s(sqrt(CEddyDist0),bs="cs",k=5),
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
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'_Daily.png',sep=""),width=600,height=600)
    plot.gam(siteDayModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
    while (dev.cur()>1) {dev.off()}
    
  }
  
  if (sum(weekData$Pres>0)>25){
    #run full weekly model for this site
    fullSiteWeekMod = gam(Pres ~ s(Sal700,bs="cs",k=5)
                          + s(Temp0,bs="cs",k=5)
                          + s(sqrt(VelAsp700),bs="cc",k=5)
                          + s(VelMag700,bs="cs",k=5)
                          + s(sqrt(AEddyDist0),bs="cs",k=5)
                          + s(sqrt(CEddyDist0),bs="cs",k=5),
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
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'_Weekly.png',sep=""),width=600,height=600)
    plot.gam(siteWeekModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
    while (dev.cur()>1) {dev.off()}
  }
}

save(siteDayModList,pValDayList,siteDayModCompList,file=paste(outDir,'/',spec,'/','DailySiteSpecificModels.Rdata',sep=""))
save(siteWeekModList,pValWeekList,siteWeekModCompList,file=paste(outDir,'/',spec,'/','WeeklySiteSpecificModels.Rdata',sep=""))
