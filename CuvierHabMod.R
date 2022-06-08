library(tidyverse)
# library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

## GAM approach ---------------------
# Regional model
spec = 'Cuvier'
outDir = "E:/ModelingCovarData/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/Cuvier_masterDF.csv'))
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
smoothVarList = c ("FSLE0",
                   "Sal0",
                   "Sal700",
                   "SSH",
                   "Temp0",
                   "Temp700",
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
# FSLE0      "140340.455536642" "135693.090137227" "134257.321481422" "134190.159811598" "fiveKnots"
# Sal0       "150936.489983644" "150863.572262892" "136325.840951809" "129891.614577069" "fiveKnots"
# Sal700     "151639.461000963" "151507.8897522"   "133712.067489804" "127212.095946889" "fiveKnots"
# SSH        "164393.149675362" "94221.3144218556" "82653.6023318413" "81913.847321494"  "fiveKnots"
# Temp0      "156592.58379896"  "156587.15460589"  "155255.399264108" "155227.206313321" "fiveKnots"
# Temp700    "152811.844010939" "147412.200061237" "135797.491181525" "136432.870583687" "fourKnots"
# CEddyDist0 "153869.762369917" "151972.080071443" "150973.847651684" "150324.10725468"  "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=4),
              data=data,
              family=poisson,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=4),
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

#                     para s(sqrt(CEddyDist0)) s(log(abs(FSLE0))) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700)
# para                   1              0.0000             0.0000  0.0000    0.0000  0.0000   0.0000     0.0000
# s(sqrt(CEddyDist0))    0              1.0000             0.0339  0.1544    0.1578  0.1012   0.0673     0.2633
# s(log(abs(FSLE0)))     0              0.0403             1.0000  0.1152    0.1119  0.1212   0.0470     0.2185
# s(Sal0)                0              0.1665             0.1061  1.0000    0.6566  0.2936   0.2464     0.8307
# s(Sal700)              0              0.1669             0.1096  0.6580    1.0000  0.3471   0.2254     0.8767
# s(SSH0)                0              0.1550             0.1005  0.3984    0.3781  1.0000   0.2985     0.7828
# s(Temp0)               0              0.1003             0.0556  0.3325    0.3239  0.2669   1.0000     0.5091
# s(Temp700)             0              0.1581             0.1068  0.5377    0.6309  0.3305   0.3347     1.0000


# Sal0 concurved with Sal700 and Temp700
# Sal700 concurved with Sal0 and Temp700
# Temp700 concurved with Sal700, Sal0, SSH, Temp0  
# check pvals of Sal0, Temp0, Sal700, Temp700, SSH

# check p vals of concurved covars
format.pval(summary(gam(Pres~s(Sal0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Sal700,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Temp700,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(Temp0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16
format.pval(summary(gam(Pres~s(SSH0,bs="cs",k=5),data=data,family=poisson,gamma=1.4,method="REML",select=TRUE))$s.pv,digits=10)
# < 2.220446e-16

# all p vlaues are the same :(
# remove temp 700 because it is concurved with sal0 & sal700 & SSH0
# remove sal0 because still still concurved w sal700


# run reduced model
dayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
             + s(log(abs(FSLE0)),bs="cs",k=5)
             # + s(Sal0,bs="cs",k=5)
             + s(Sal700,bs="cs",k=5)
             + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5),
             # + s(Temp700,bs="cs",k=4),
             data=data,
             family=poisson,
             method="REML",
             select=TRUE,
             gamma=1.4,
             na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              # + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5),
              # + s(Temp700,bs="cs",k=4),
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

#                     para s(sqrt(CEddyDist0)) s(log(abs(FSLE0))) s(Sal700) s(SSH0) s(Temp0)
# para                   1              0.0000             0.0000    0.0000  0.0000   0.0000
# s(sqrt(CEddyDist0))    0              1.0000             0.0554    0.1902  0.1092   0.0814
# s(log(abs(FSLE0)))     0              0.0671             1.0000    0.1727  0.1837   0.0653
# s(Sal700)              0              0.2016             0.1782    1.0000  0.4132   0.2330
# s(SSH0)                0              0.1840             0.1607    0.4023  1.0000   0.3026
# s(Temp0)               0              0.1216             0.0866    0.3179  0.2631   1.0000

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
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(Sal700, 
#                                                   bs = "cs", k = 5) + s(sqrt(CEddyDist0), bs = "cs", 
#                                                                         k = 5) + s(SSH0, bs = "cs", k = 5) + s(Temp0, bs = "cs", 
#                                                                                                                k = 5) + 1
# 
# Parametric coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.69079    0.03005  -22.99   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                           edf Ref.df  Chi.sq p-value    
#   s(log(abs(FSLE0)))  3.982      4   796.3  <2e-16 ***
#   s(Sal700)           3.990      4  7204.9  <2e-16 ***
#   s(sqrt(CEddyDist0)) 3.550      4   229.6  <2e-16 ***
#   s(SSH0)             3.994      4 12861.7  <2e-16 ***
#   s(Temp0)            3.957      4   836.6  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =   0.47   Deviance explained =   63%
# -REML =  24358  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(Sal700, 
#                                                   bs = "cs", k = 5) + s(sqrt(CEddyDist0), bs = "cs", 
#                                                                         k = 5) + s(SSH0, bs = "cs", k = 5) + s(Temp0, bs = "cs", 
#                                                                                                                k = 5) + 1
# 
# Parametric coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.25379    0.03046   41.17   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                           edf Ref.df  Chi.sq p-value    
#   s(log(abs(FSLE0)))  3.971      4   367.9  <2e-16 ***
#   s(Sal700)           3.993      4  6407.2  <2e-16 ***
#   s(sqrt(CEddyDist0)) 3.956      4   656.2  <2e-16 ***
#   s(SSH0)             3.993      4 12131.6  <2e-16 ***
#   s(Temp0)            3.961      4  1020.4  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.638   Deviance explained = 74.9%
# -REML =  13780  Scale est. = 1         n = 1509

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
                         + s(log(abs(FSLE0)),bs="cs",k=5)
                         # + s(Sal0,bs="cs",k=5)
                         + s(Sal700,bs="cs",k=5)
                         + s(SSH0,bs="cs",k=5)
                         + s(Temp0,bs="cs",k=5),
                         # + s(Temp700,bs="cs",k=4),
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
                          + s(log(abs(FSLE0)),bs="cs",k=5)
                          # + s(Sal0,bs="cs",k=5)
                          + s(Sal700,bs="cs",k=5)
                          + s(SSH0,bs="cs",k=5)
                          + s(Temp0,bs="cs",k=5),
                          # + s(Temp700,bs="cs",k=4),
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

