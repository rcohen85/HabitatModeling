library(tidyverse)
library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

# library(tidyverse)
# library(splines2)
# library(geepack)
# source("getPvalues.R")

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/Sowerby_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

## GAM approach ---------------------
# Regional model
spec = 'Sowerby'
outDir = "E:/ModelingCovarData/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/Sowerby_masterDF.csv'))
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
smoothVarList = c("Chl0",
                  "FSLE0",
                  "Sal0",
                  "Sal700",
                  "SSH0",
                  "Temp0",
                  "Temp700",
                  "VelMag0")


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

#             linMod             threeKnots         fourKnots          fiveKnots          Best       
# Chl0    "35013.0853510471" "33232.0232032764" "31660.6009127515" "31627.5984841782" "fiveKnots"
# FSLE0   "35702.6922112411" "35644.7729250646" "35520.7095082101" "35512.1108657343" "fiveKnots"
# Sal0    "32622.5654441784" "30482.403230816"  "29142.9869296782" "29147.341221908"  "fourKnots"
# Sal700  "32629.4093685009" "30861.0270009168" "29003.3655042682" "28939.4633678373" "fiveKnots"
# SSH0    "29231.9517730921" "28840.4512187974" "28833.7218307997" "28815.1116948551" "fiveKnots"
# Temp0   "32829.2851872648" "32387.1979406793" "32298.3754511516" "32240.0186657552" "fiveKnots"
# Temp700 "30162.5483503848" "29150.0275790315" "28911.9937925605" "28909.898174784"  "fiveKnots"
# VelMag0 "35046.5738366041" "34748.686339373"  "34741.8234827317" "34628.5883860314" "fiveKnots"


# run full model
fullMod = gam(Pres ~ s(VelMag0,bs="cs",k=5)
              + s(log(Chl0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=5),
              data=data,
              family=poisson,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

weekMod = gam(Pres ~ s(VelMag0,bs="cs",k=5)
              + s(log(Chl0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=5),
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

#                     para s(VelMag0) s(log(Chl0)) s(log(abs(FSLE0))) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700)
# para                  1     0.0000       0.0000             0.0000  0.0000    0.0000  0.0000   0.0000     0.0000
# s(VelMag0)            0     1.0000       0.0322             0.0522  0.2141    0.1537  0.2403   0.0734     0.1375
# s(log(Chl0))          0     0.0391       1.0000             0.0429  0.2874    0.2208  0.2121   0.2897     0.2419
# s(log(abs(FSLE0)))    0     0.1094       0.0481             1.0000  0.1525    0.1119  0.1212   0.0470     0.1250
# s(Sal0)               0     0.1799       0.2629             0.0938  1.0000    0.6222  0.2479   0.2140     0.5613
# s(Sal700)             0     0.2893       0.2947             0.1096  0.7501    1.0000  0.3471   0.2254     0.6630
# s(SSH0)               0     0.3099       0.3703             0.1005  0.4963    0.3781  1.0000   0.2985     0.4560
# s(Temp0)              0     0.1362       0.3724             0.0556  0.3451    0.3239  0.2669   1.0000     0.3935
# s(Temp700)            0     0.1966       0.2989             0.1076  0.6538    0.6750  0.3340   0.3463     1.0000

# VelMag0 concurved with Sal700, SSH
# Chl0 concurved with Sal0, Sal700, SSH, Temp0, Temp700
# Sal0 concurved with Chl0, Sal700, SSH0, Temp0, Temp700
# Sal700 concurved with Sal0, SSH, Temp0, and Temp700
# Temp700 concurved with Sal700
# SSH0 concurved with Sal700, Temp0, Temp700
# Temp0 concurved with Chl0, SSH, Temp700
# Temp700 concurved with Sal0, Sal700, SSH, Temp0
# check pvals of VelMag0, Chl0, SSH, Sal0, Sal700, Temp0, Temp700

# check p vals of concurved covars
format.pval(summary(gam(Pres~s(VelMag0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Chl0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(SSH0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Sal0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Sal700,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Temp0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Temp700,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16

# all p vlaues are the same :(
# remove SSH because it is concurved with literally everything
# remove temp700, sal0 because concurved with everything


# run reduced model
dayMod = gam(Pres ~ s(VelMag0,bs="cs",k=5)
             + s(log(Chl0),bs="cs",k=5)
             + s(log(abs(FSLE0)),bs="cs",k=5)
             # + s(Sal0,bs="cs",k=4)
             + s(Sal700,bs="cs",k=5)
             # + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5),
             # + s(Temp700,bs="cs",k=5),
             data=data,
             family=poisson,
             method="REML",
             select=TRUE,
             gamma=1.4,
             na.action="na.fail")

weekMod = gam(Pres ~ s(VelMag0,bs="cs",k=5)
              + s(log(Chl0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              # + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              # + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5),
              # + s(Temp700,bs="cs",k=5),
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

#                     para s(VelMag0) s(log(Chl0)) s(log(abs(FSLE0))) s(Sal700) s(Temp0)
# para                  1     0.0000       0.0000             0.0000    0.0000   0.0000
# s(VelMag0)            0     1.0000       0.0637             0.0882    0.2204   0.0994
# s(log(Chl0))          0     0.0758       1.0000             0.0723    0.2177   0.2900
# s(log(abs(FSLE0)))    0     0.1353       0.0759             1.0000    0.1727   0.0653
# s(Sal700)             0     0.3869       0.3001             0.1782    1.0000   0.2330
# s(Temp0)              0     0.1811       0.3618             0.0866    0.3179   1.0000

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
#   Pres ~ s(log(Chl0), bs = "cs", k = 5) + s(Sal700, bs = "cs", 
#                                             k = 5) + s(Temp0, bs = "cs", k = 5) + s(VelMag0, bs = "cs", 
#                                                                                     k = 5) + 1
# 
# Parametric coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.93829    0.07425  -26.11   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                   edf Ref.df Chi.sq p-value    
#   s(log(Chl0)) 2.626      4 178.72  <2e-16 ***
#   s(Sal700)    3.982      4 667.39  <2e-16 ***
#   s(Temp0)     3.886      4 196.56  <2e-16 ***
#   s(VelMag0)   3.581      4  40.31  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.112   Deviance explained = 28.7%
# -REML =  10187  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(log(Chl0), 
#                                                   bs = "cs", k = 5) + s(Sal700, bs = "cs", k = 5) + 
#   s(Temp0, bs = "cs", k = 5) + s(VelMag0, bs = "cs", 
#                                  k = 5) + 1
# 
# Parametric coefficients:
#              Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.2660     0.1012  -2.629  0.00857 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                         edf Ref.df Chi.sq p-value    
#   s(log(abs(FSLE0))) 3.861      4 163.47  <2e-16 ***
#   s(log(Chl0))       2.772      4 128.38  <2e-16 ***
#   s(Sal700)          3.976      4 507.56  <2e-16 ***
#   s(Temp0)           3.866      4 198.96  <2e-16 ***
#   s(VelMag0)         3.659      4  52.38  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.201   Deviance explained = 42.1%
# -REML = 5477.9  Scale est. = 1         n = 1509
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
    fullSiteDayMod = gam(Pres ~ s(VelMag0,bs="cs",k=5)
                         + s(log(Chl0),bs="cs",k=5)
                         + s(log(abs(FSLE0)),bs="cs",k=5)
                         # + s(Sal0,bs="cs",k=4)
                         + s(Sal700,bs="cs",k=5)
                         # + s(SSH0,bs="cs",k=5)
                         + s(Temp0,bs="cs",k=5),
                         # + s(Temp700,bs="cs",k=5),
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
    fullSiteWeekMod = gam(Pres ~ s(VelMag0,bs="cs",k=5)
                          + s(log(Chl0),bs="cs",k=5)
                          + s(log(abs(FSLE0)),bs="cs",k=5)
                          # + s(Sal0,bs="cs",k=4)
                          + s(Sal700,bs="cs",k=5)
                          # + s(SSH0,bs="cs",k=5)
                          + s(Temp0,bs="cs",k=5),
                          # + s(Temp700,bs="cs",k=5),
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

