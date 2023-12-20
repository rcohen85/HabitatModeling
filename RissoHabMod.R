library(tidyverse)
library(stringr)
library(mgcv)
library(splines)
library(AER)
library(MuMIn)
library(pracma)

## GAM approach ---------------------
# Regional model
spec = 'Risso'
outDir = "E:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('Risso_masterDF_plusJAX.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)
data$Date = as.Date(data$Date,origin='1970-01-01')

# create weekly time series to reduce autocorrelation
stDt = as.Date("2016-05-01")
edDt = as.Date("2019-04-30")
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
allDates = seq.Date(stDt,edDt,by=1)
weekDates = seq.Date(stDt,edDt,by=7)
weeklyDF = as.numeric()

for (l in 1:length(sites)) {
  
  # Create dataframe to hold data (or NAs) for all dates
  fullDatesDF = data.frame(matrix(nrow=length(allDates), ncol=dim(data)[2]))
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
  
  # Sum presence in each week
  summaryData = fullDatesDF %>%
    group_by(WeekID) %>%
    summarize(Pres=sum(Pres,na.rm=TRUE))
  
  # normalize by effort
  effDF = data.frame(count=rep(1,length(allDates)),weekID=weekID)
  effDF$count[which(is.na(fullDatesDF$Pres))] = 0
  propEff = effDF %>%
    group_by(weekID) %>%
    summarize(sum(count))
  summaryData$propEff = unlist(propEff[,2])/7
  summaryData$Pres = summaryData$Pres*(1/summaryData$propEff)
  
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

# Center and scale all predictors
weeklyDF[,6:19] = scale(weeklyDF[,6:19],center=TRUE,scale=TRUE)

# Remove incomplete observations (NAs in FSLE)
badRows = unique(which(is.na(weeklyDF),arr.ind=TRUE)[,1])
weeklyDF = weeklyDF[-badRows,]

# re-round presence data
weeklyDF$Pres = round(weeklyDF$Pres)

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("CEddyDist0",
                  "Chl0",
                  "EKE0",
                  "FSLE0",
                  "Sal0",
                  "SSH0",
                  "Temp0")

# check residual autocorrelation of weekly data
sites = unique(weeklyDF$Site)
residAutocorr = matrix(ncol=1,nrow=length(sites))
rownames(residAutocorr) = sites
for (j in 1:length(sites)){
  
  siteInd = which(!is.na(str_match(weeklyDF$Site,sites[j])))
  
  if (sum(which(weeklyDF$Pres[siteInd]>0))>10){
    
    siteData = weeklyDF[siteInd,]
    
    BlockMod = glm(Pres~ bs(CEddyDist0) # include all terms in smoothVarList above!!
                   + bs(Chl0)
                   + bs(EKE0)
                   + bs(FSLE0)
                   + bs(Sal0)
                   + bs(SSH0)
                   + bs(Temp0),
                   data=siteData,family=poisson)
    
    acorr = acf(residuals(BlockMod), lag.max = 1000, main=paste(spec,"at",sites[j]))
    CI = ggfortify:::confint.acf(acorr)
    ACFidx = which(acorr[["acf"]] < CI, arr.ind=TRUE)
    residAutocorr[j,1] = ACFidx[1]
  }
}
residAutocorr

# HZ     5
# OC     5
# NC     4
# BC     3
# WC     3
# NFC    8
# HAT    2
# GS     3
# BP     3
# BS     2
# JAX    5

# test for overdispersion
dispMod = glm(Pres~ bs(CEddyDist0) # include all terms in smoothVarList above!!
              + bs(Chl0)
              + bs(EKE0)
              + bs(FSLE0)
              + bs(Sal0)
              + bs(SSH0)
              + bs(Temp0),
              data=weeklyDF,family=poisson)

dispersiontest(dispMod,alternative='two.sided')
# Dispersion test
# 
# data:  dispMod
# z = 8.171, p-value = 3.058e-16
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 92.58225 

# data are very overdispersed, will use Tweedie family in models
modFam=tw

# run full model
weekMod = gam(Pres ~ + s(CEddyDist0,bs="ts",k=4)
              + s(Chl0,bs="ts",k=4)
              + s(EKE0,bs="ts",k=4)
              + s(FSLE0,bs="ts",k=4)
              + s(Sal0,bs="ts",k=4)
              + s(SSH0,bs="ts",k=4)
              + s(Temp0,bs="ts",k=4),
              data=weeklyDF,
              family=modFam,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              optimizer=c('outer','bfgs'),
              select=TRUE)

# check convergence
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                 para s(CEddyDist0) s(Chl0) s(EKE0) s(FSLE0) s(Sal0) s(SSH0) s(Temp0)
# para             1        0.0000  0.0000  0.0000   0.0000  0.0000  0.0000   0.0000
# s(CEddyDist0)    0        1.0000  0.0404  0.0230   0.0346  0.1824  0.1974   0.1164
# s(Chl0)          0        0.1138  1.0000  0.0113   0.0259  0.3332  0.6354   0.5345
# s(EKE0)          0        0.0403  0.0086  1.0000   0.0366  0.1266  0.0496   0.0377
# s(FSLE0)         0        0.0304  0.0457  0.0259   1.0000  0.2263  0.1687   0.1374
# s(Sal0)          0        0.1899  0.1790  0.0449   0.1182  1.0000  0.6840   0.4449
# s(SSH0)          0        0.1906  0.2239  0.0298   0.1005  0.4942  1.0000   0.6452
# s(Temp0)         0        0.0989  0.2647  0.0211   0.0694  0.3406  0.4541   1.0000
# > 

# SSH0 concurved w Sal0, Chl0
# Temp0 concurved w Chl0, SSH0

# taking out SSH0

weekMod = gam(Pres ~ + s(CEddyDist0,bs="ts",k=4)
              + s(Chl0,bs="ts",k=4)
              + s(EKE0,bs="ts",k=4)
              + s(FSLE0,bs="ts",k=4)
              + s(Sal0,bs="ts",k=4)
              # + s(SSH0,bs="ts",k=4)
              + s(Temp0,bs="ts",k=4),
              data=weeklyDF,
              family=modFam,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              optimizer=c('outer','bfgs'),
              select=TRUE)

conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                 para s(CEddyDist0) s(Chl0) s(EKE0) s(FSLE0) s(Sal0) s(Temp0)
# para             1        0.0000  0.0000  0.0000   0.0000  0.0000   0.0000
# s(CEddyDist0)    0        1.0000  0.0404  0.0230   0.0346  0.1824   0.1164
# s(Chl0)          0        0.1138  1.0000  0.0113   0.0259  0.3332   0.5345
# s(EKE0)          0        0.0403  0.0086  1.0000   0.0366  0.1266   0.0377
# s(FSLE0)         0        0.0304  0.0457  0.0259   1.0000  0.2263   0.1374
# s(Sal0)          0        0.1899  0.1790  0.0449   0.1182  1.0000   0.4449
# s(Temp0)         0        0.0989  0.2647  0.0211   0.0694  0.3406   1.0000

# Temp0 still convurved w Chl0
# taking out Chl0

weekMod = gam(Pres ~ + s(CEddyDist0,bs="ts",k=4)
              # + s(Chl0,bs="ts",k=4)
              + s(EKE0,bs="ts",k=4)
              + s(FSLE0,bs="ts",k=4)
              + s(Sal0,bs="ts",k=4)
              # + s(SSH0,bs="ts",k=4)
              + s(Temp0,bs="ts",k=4),
              data=weeklyDF,
              family=modFam,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              optimizer=c('outer','bfgs'),
              select=TRUE)

conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                 para s(CEddyDist0) s(EKE0) s(FSLE0) s(Sal0) s(Temp0)
# para             1        0.0000  0.0000   0.0000  0.0000   0.0000
# s(CEddyDist0)    0        1.0000  0.0230   0.0346  0.1824   0.1164
# s(EKE0)          0        0.0403  1.0000   0.0366  0.1266   0.0377
# s(FSLE0)         0        0.0304  0.0259   1.0000  0.2263   0.1374
# s(Sal0)          0        0.1899  0.0449   0.1182  1.0000   0.4449
# s(Temp0)         0        0.0989  0.0211   0.0694  0.3406   1.0000

# no convurvity >0.5

# dayModCompTable = dredge(dayMod,
#                          beta="none",
#                          evaluate=TRUE,
#                          trace=TRUE)

weekModCompTable = dredge(weekMod,
                          beta="none",
                          rank='AIC',
                          evaluate=TRUE,
                          trace=TRUE)

# Save best models
# optWeekMod = get.models(weekModCompTable,subset=1)
# optWeekMod = optWeekMod[[names(optWeekMod)]]
# save(optWeekMod,weekModCompTable,file=paste(outDir,'/',spec,'/','WeeklyRegionalModel.Rdata',sep=""))
topMods = get.models(weekModCompTable,subset=delta<2) # all best-performing models
if (numel(names(topMods))>1){
  optWeekMod = model.avg(weekModCompTable,subset=delta<2,fit=TRUE) # averaged model
  save(optWeekMod,topMods,weekModCompTable,file=paste(spec,'_','WeeklyRegionalModel_Updated.Rdata',sep=""))
}else{optWeekMod = topMods[[1]]
save(optWeekMod,weekModCompTable,file=paste(spec,'_','WeeklyRegionalModel_Updated.Rdata',sep=""))}

# sink(paste(outDir,'/',spec,'/','WeeklyRegionalModelSummary.txt',sep=""))
sink(paste(spec,'_','WeeklyRegionalModelSummary_Updated.txt',sep=""))
for (i in 1:length(names(topMods))){
  print(summary(topMods[[names(topMods)[i]]]))}
sink()

# plot
# png(filename=paste(spec,'_allSitesWeekly_Updated.png',sep=""),width=600,height=600)
# plot.gam(optWeekMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0,main="Optimal Model")
# while (dev.cur()>1) {dev.off()}

# Retrain model with 2/3 of data, then validate by predicting remaining 1/3
trainInd = sample(1:nrow(weeklyDF),floor(nrow(weeklyDF)*.66))
testInd = setdiff(1:nrow(weeklyDF),trainInd)

if (numel(names(topMods))>1){
  valModList = list()
  for (i in 1:numel(names(topMods))){
    valMod = update(topMods[[1]],data=weeklyDF[trainInd,])
    valModList[[i]] = valMod
  }
  valMod = model.avg(valModList,subset=delta<2,fit=TRUE)
  save(valModList,valMod,file=paste(spec,'_','ValidationModel_Updated.Rdata',sep=""))
}else{valMod = update(optWeekMod,data=weeklyDF[trainInd,])
save(valMod,file=paste(spec,'_','ValidationModel_Updated.Rdata',sep=""))}

true = weeklyDF[testInd,"Pres"]$Pres
preds = predict(valMod,weeklyDF[testInd,],full=TRUE,type="response",backtransform=FALSE)
plot(true,preds,type="p")
cor(true,preds)


# # Site-specific models ---------------
# sites = unique(data$Site)
# siteDayModList = list()
# pValDayList = list()
# siteDayModCompList = list()
# siteWeekModList = list()
# pValWeekList = list()
# siteWeekModCompList = list()
# for (i in 1:length(sites)){
#   
#   dayInd = which(!is.na(str_match(data$Site,sites[i])))
#   dayData = data[dayInd,]
#   weekInd = which(!is.na(str_match(weeklyDF$Site,sites[i])))
#   weekData = weeklyDF[weekInd,]
#   
#   if (sum(dayData$Pres>0)>25){
#     
#     # run full daily model for this site
#     fullSiteDayMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
#                          + s(log(abs(FSLE0)),bs="cs",k=5)
#                          + s(Sal400,bs="cs",k=5)
#                          + s(Temp0,bs="cs",k=5)
#                          + s(VelMag400,bs="cc",k=5)
#                          + s(sqrt(AEddyDist0),bs="cs",k=5)
#                          + s(sqrt(CEddyDist0),bs="cs",k=5),
#                          data=dayData,
#                          family=poisson,
#                          method="REML",
#                          select=TRUE,
#                          gamma=1.4,
#                          na.action="na.fail")
#     
#     siteDayModCompTable = dredge(fullSiteDayMod,
#                                  beta="none",
#                                  evaluate=TRUE,
#                                  trace=TRUE)
#     siteDayModCompList[[sites[i]]] = siteDayModCompTable
#     
#     optSiteDayMod = get.models(siteDayModCompTable,subset=1)
#     optSiteDayMod = optSiteDayMod[[names(optSiteDayMod)]]
#     siteDayPV = summary(optSiteDayMod)$s.pv
#     
#     if (any(siteDayPV>=0.05)){ # Remove non-significant terms & re-run model iteratively until only signif covars remain
#       flag = 1
#       while (flag==1){
#         # get terms from formula as strings
#         thisForm = as.character(optSiteDayMod$formula)[3]
#         startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
#         termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
#         termInd = c(0,termInd,str_length(thisForm)+1)
#         allTerms = character()
#         for (j in 1:length(termInd)-1){
#           thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
#           allTerms = c(allTerms,thisTerm)
#         }
#         # identify which terms were non-significant
#         badVars = allTerms[siteDayPV>=0.05]
#         dontNeed = which(!is.na(str_match(badVars,"1")))
#         if (!is_empty(dontNeed)){
#           badVars = badVars[-dontNeed]}
#         # update model
#         optSiteDayMod<-eval(parse(text=paste("update(optSiteDayMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
#         siteDayPV = summary(optSiteDayMod)$s.pv
#         if (!any(siteDayPV>=0.05)){
#           siteDayModList[[sites[i]]] = optSiteDayMod
#           pValDayList[[sites[i]]] = siteDayPV
#           flag=0
#         }
#       }
#     } else {
#       siteDayModList[[sites[i]]] = optSiteDayMod
#       pValDayList[[sites[i]]] = siteDayPV
#     }
#     
#     sink(paste(outDir,'/',spec,'/',spec,'_',sites[i],'_DailySummary.txt',sep=""))
#     print(summary(siteDayModList[[sites[i]]]))
#     sink()
#     
#     png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'_Daily.png',sep=""),width=600,height=600)
#     plot.gam(siteDayModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
#     while (dev.cur()>1) {dev.off()}
#     
#   }
#   
#   if (sum(weekData$Pres>0)>25){
#     #run full weekly model for this site
#     fullSiteWeekMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
#                           + s(log(abs(FSLE0)),bs="cs",k=5)
#                           + s(Sal400,bs="cs",k=5)
#                           + s(Temp0,bs="cs",k=5)
#                           + s(VelMag400,bs="cc",k=5)
#                           + s(sqrt(AEddyDist0),bs="cs",k=5)
#                           + s(sqrt(CEddyDist0),bs="cs",k=5),
#                           data=weekData,
#                           family=poisson,
#                           method="REML",
#                           select=TRUE,
#                           gamma=1.4,
#                           na.action="na.fail")
#     
#     siteWeekModCompTable = dredge(fullSiteWeekMod,
#                                   beta="none",
#                                   evaluate=TRUE,
#                                   trace=TRUE)
#     siteWeekModCompList[[sites[i]]] = siteWeekModCompTable
#     
#     optSiteWeekMod = get.models(siteWeekModCompTable,subset=1)
#     optSiteWeekMod = optSiteWeekMod[[names(optSiteWeekMod)]]
#     siteWeekPV = summary(optSiteWeekMod)$s.pv
#     
#     if (any(siteWeekPV>=0.05)){ # Remove non-significant terms & re-run model iteratively until only signif covars remain
#       flag = 1
#       while (flag==1){
#         # get terms from formula as strings
#         thisForm = as.character(optSiteWeekMod$formula)[3]
#         startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
#         termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
#         termInd = c(0,termInd,str_length(thisForm)+1)
#         allTerms = character()
#         for (j in 1:length(termInd)-1){
#           thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
#           allTerms = c(allTerms,thisTerm)
#         }
#         # identify which terms were non-significant
#         badVars = allTerms[siteWeekPV>=0.05]
#         dontNeed = which(!is.na(str_match(badVars,"1")))
#         if (!is_empty(dontNeed)){
#           badVars = badVars[-dontNeed]}
#         # update model
#         optSiteWeekMod<-eval(parse(text=paste("update(optSiteWeekMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
#         siteWeekPV = summary(optSiteWeekMod)$s.pv
#         if (!any(siteWeekPV>=0.05)){
#           siteWeekModList[[sites[i]]] = optSiteWeekMod
#           pValWeekList[[sites[i]]] = siteWeekPV
#           flag=0
#         }
#       }
#     } else {
#       siteWeekModList[[sites[i]]] = optSiteWeekMod
#       pValWeekList[[sites[i]]] = siteWeekPV
#     }
#     
#     sink(paste(outDir,'/',spec,'/',spec,'_',sites[i],'_WeeklySummary.txt',sep=""))
#     print(summary(siteWeekModList[[sites[i]]]))
#     sink()
#     
#     png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'_Weekly.png',sep=""),width=600,height=600)
#     plot.gam(siteWeekModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
#     while (dev.cur()>1) {dev.off()}
#   }
# }
# 
# save(siteDayModList,pValDayList,siteDayModCompList,file=paste(outDir,'/',spec,'/','DailySiteSpecificModels.Rdata',sep=""))
# save(siteWeekModList,pValWeekList,siteWeekModCompList,file=paste(outDir,'/',spec,'/','WeeklySiteSpecificModels.Rdata',sep=""))
# 
# 
# 
# 
# 
# # GEEGLM approach ------------------------------------
# lagID = 43
# numClust = length(data$Pres)/(lagID-1)
# if (numClust<length(data$Pres)){
#   clustID = rep(1:ceiling(numClust),each=lagID)
#   clustID = clustID[1:length(data$Pres)]
# } else {
#   clustID = 1:length(data$Pres)
# }
# data$GroupID = clustID
# 
# # Round presence to get Poisson dist
# data$Pres = round(data$Pres)
# 
# # Test for how a term should be included in the model
# startTime = Sys.time()
# smoothVarList = c("Temp0",
#                   "Sal0",
#                   "Sal200",
#                   "SSH0",
#                   "FSLE0",
#                   "GSLat",
#                   "Slope",
#                   "Aspect")
# 
# modOpts = c("linMod","threeKnots","fourKnots")
# QIC_votes = matrix(nrow=length(smoothVarList),ncol=4)
# 
# for (i in 1:(length(smoothVarList)-1)){
#   
#   modelCall = paste('geeglm(Pres~data$',smoothVarList[i],',data=data,family=poisson,id=GroupID,corstr="ar1")',sep="")
#   linMod = eval(parse(text=modelCall))
#   
#   modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.5)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],'))),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
#   smoothMod1 = eval(parse(text=modelCall))
#   
#   modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.333,0.666)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],'))),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
#   smoothMod2 = eval(parse(text=modelCall))
#   
#   # modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.275,0.5,0.725)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],'))),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
#   # smoothMod3 = eval(parse(text=modelCall))
#   
#   QIC_votes[i,1:3] = c(QIC(linMod)[[1]],QIC(smoothMod1)[[1]],QIC(smoothMod2)[[1]])
#   QIC_votes[i,4] = modOpts[which.min(QIC_votes[i,1:3])]
# }
# 
# endTime = Sys.time()
# endTime-startTime
# colnames(QIC_votes) = c(modOpts,"Best")
# rownames(QIC_votes) = smoothVarList[]
# QIC_votes
# 
# #           linMod              threeKnots          fourKnots           Best        
# # Temp0  "-2458133.44127198" "-2497958.76610649" "-2494447.98398555" "threeKnots"
# # Sal0   "-2387145.82259727" "-2397529.44160543" "-2397770.58164437" "fourKnots" 
# # Sal200 "-2385452.46919174" "-2396045.24119207" "-2397783.78049046" "fourKnots" 
# # SSH0   "-2509932.72715967" "-2510614.47948344" "-2512094.69825831" "fourKnots" 
# # FSLE0  "-2385328.82959341" "-2386058.84470483" "-2386081.53007312" "fourKnots" 
# # GSLat  "-2386298.16752796" "-2387378.29830821" "-2387372.85097507" "threeKnots"
# # Slope  "-2404596.46021069" "-2530688.29829264" "-2532163.73313361" "fourKnots" 
# # Aspect NA                  NA                  NA                  NA          
# 
# # Make smooth terms, run full model and check collinearity
# smoothVarList = c("Temp0",
#                   "Sal0",
#                   "Sal200",
#                   "SSH0",
#                   "FSLE0",
#                   "GSLat",
#                   "Slope",
#                   "Aspect")
# knotList = list(c(0.5),
#                 c(0.333,0.666),
#                 c(0.333,0.666),
#                 c(0.333,0.666),
#                 c(0.333,0.666),
#                 c(0.5),
#                 c(0.333,0.666),
#                 c(0.333,0.666))
# linVarList = list()
# smoothNameList = character()
# 
# 
# for (i in 1:length(smoothVarList)){
#   
#   if (str_detect(smoothVarList[i],"Asp")){
#     eval(parse(text=paste('S_',smoothVarList[i],'= mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=unlist(knotList[i])),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')),periodic=TRUE)',sep="")))
#   } else {
#     eval(parse(text=paste('S_',smoothVarList[i],'= mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=unlist(knotList[i])),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')))',sep="")))
#   }
#   
#   smoothNameList = c(smoothNameList,paste('S_',smoothVarList[i],sep=""))
# }
# thisForm = formula(paste('Pres~',paste(c(smoothNameList,linVarList),collapse="+"),sep=""))
# 
# fullMod = geeglm(thisForm,
#                  family=poisson,
#                  data=data,
#                  id=GroupID,
#                  corstr="unstructured")
# 
# VIFvals = vif(fullMod)
# VIFvals = cbind(VIFvals,(VIFvals[,3])^2)
# colnames(VIFvals)[4] = "LOOK AT ME"
# VIFvals
# 
# #               GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# # S_Temp0   8.585876  4        1.308347   1.711773
# # S_Sal0   31.063539  5        1.410019   1.988155
# # S_Sal200 34.562502  5        1.425150   2.031052
# # S_SSH0   49.168005  5        1.476278   2.179398
# # S_FSLE0   4.956267  5        1.173587   1.377308
# # S_GSLat   2.698404  4        1.132109   1.281672
# # S_Slope  88.670559  5        1.565950   2.452200
# # S_Aspect 11.408530  2        1.837839   3.377652
# 
# # no collinearity, no need to remove any terms
# 
# # check convergence
# fullMod$geese$error
# # 1 model didn't converge
# 
# # run w. unstructured correlation
# fullMod = geeglm(thisForm,
#                  family=poisson,
#                  data=data,
#                  id=GroupID,
#                  corstr="unstructured")
# fullMod$geese$error
# # 0 model converged
# 
# #               GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# # S_Temp0   5.896675  4        1.248320   1.558303
# # S_Sal0   14.987336  5        1.310909   1.718482
# # S_Sal200 16.682475  5        1.325031   1.755707
# # S_SSH0   30.147910  5        1.405807   1.976293
# # S_FSLE0   3.262515  5        1.125525   1.266807
# # S_GSLat   2.855038  4        1.140122   1.299879
# # S_Slope  78.268675  5        1.546532   2.391760
# # S_Aspect 12.849682  2        1.893316   3.584645
# 
# # check term significance
# PV = getPvalues(fullMod)
# 
# #Remove non-significant terms, re-run model, check p-values
# PV$'p-value'[PV$'p-value'=="<0.0001"] = 0.0001
# badVars = PV$Variable[as.numeric(PV$'p-value')>=0.05]
# redMod<-eval(parse(text=paste("update(fullMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
# if (redMod$geese$error==1){
#   print("Model did not converge")
# } else {PVred = getPvalues(redMod)
# PVred$'p-value'[PVred$'p-value'=="<0.0001"] = 0.0001
# PVred}
# 
# # Get p-values
# PV = getPvalues(reducedMod)
# 
# # Plot terms from regional model
# source("plotSmooths.R")
# source("plotLinears.R")
# terms = names(reducedMod$model)[2:length(names(reducedMod$model))]
# k=3
# 
# for (i in 1:length(terms)){
#   if (str_detect(terms[i],"S_")){ # plot smooth terms
#     term = str_remove(terms[i],"S_")
#     coefInd = which(str_detect(names(reducedMod$coefficients),term))
#     if (str_detect(term,"Asp")){periodic=TRUE} else {periodic=FALSE}
#     print(plotSmooths(reducedMod,term,coefInd,k,periodic,site=NA,title=NULL))
#   } else { # plot linear terms
#     term=terms[i]
#     coefInd = which(str_detect(names(reducedMod$coefficients),term))
#     print(plotLinears(reducedMod,term,coefInd,site=NA,title=NULL))
#   }
# }
