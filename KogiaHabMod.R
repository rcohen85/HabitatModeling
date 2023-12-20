library(tidyverse)
library(stringr)
library(mgcv)
library(splines)
library(AER)
library(MuMIn)
library(pracma)

## GAM approach ---------------------
# Regional model
spec = 'Kogia'
outDir = "E:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('Kogia_masterDF_plusJAX.csv'))
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
  
  # sum presence in each week
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

# # don't use observations from northern sites, could be harbor porpoise
# southSites = weeklyDF$Site%in%c("HAT","GS","BP","BS")
# weeklyDF = weeklyDF[southSites,]

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c ("AEddyDist0",
                   "CEddyDist0",
                   "Chl0",
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
    
    BlockMod = glm(Pres~bs(AEddyDist0) # include all terms in smoothVarList above!!
                   + bs(CEddyDist0)
                   + bs(Chl0) 
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
# OC    NA
# NC     4
# BC     3
# WC     2
# NFC    2
# HAT    2
# GS     2
# BP     3
# BS     2
# JAX    3


# test for overdispersion
dispMod = glm(Pres~bs(AEddyDist0) # include all terms in smoothVarList above!!
              + bs(CEddyDist0)
              + bs(Chl0) 
              + bs(Sal0)
              + bs(SSH0)
              + bs(Temp0),
              data=siteData,family=poisson)

dispersiontest(dispMod,alternative='two.sided')
# Dispersion test
# 
# data:  dispMod
# z = 2.7921, p-value = 0.005237
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 4.083142     

# data are overdispersed, will use Tweedie family in models
modFam=tw

# run full model
weekMod = gam(Pres ~ s(AEddyDist0,bs="ts",k=4)
              + s(CEddyDist0,bs="ts",k=4)
              + s(Chl0,bs="ts",k=4)
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

#                 para s(AEddyDist0) s(CEddyDist0) s(Chl0) s(Sal0) s(SSH0) s(Temp0)
# para             1        0.0000        0.0000  0.0000  0.0000  0.0000   0.0000
# s(AEddyDist0)    0        1.0000        0.0058  0.0234  0.0606  0.1897   0.0234
# s(CEddyDist0)    0        0.0244        1.0000  0.0404  0.1824  0.1974   0.1164
# s(Chl0)          0        0.1060        0.1138  1.0000  0.3332  0.6354   0.5345
# s(Sal0)          0        0.1933        0.1899  0.1790  1.0000  0.6840   0.4449
# s(SSH0)          0        0.1805        0.1906  0.2239  0.4942  1.0000   0.6452
# s(Temp0)         0        0.1143        0.0989  0.2647  0.3406  0.4541   1.0000

# SSH0 concurved w Chl0, Sal0
# Temp0 concurved w Chl0, SSH0

# removing SSH0

weekMod = gam(Pres ~ s(AEddyDist0,bs="ts",k=4)
              + s(CEddyDist0,bs="ts",k=4)
              + s(Chl0,bs="ts",k=4)
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

# check convergence
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                 para s(AEddyDist0) s(CEddyDist0) s(Chl0) s(Sal0) s(Temp0)
# para             1        0.0000        0.0000  0.0000  0.0000   0.0000
# s(AEddyDist0)    0        1.0000        0.0058  0.0234  0.0606   0.0234
# s(CEddyDist0)    0        0.0244        1.0000  0.0404  0.1824   0.1164
# s(Chl0)          0        0.1060        0.1138  1.0000  0.3332   0.5345
# s(Sal0)          0        0.1933        0.1899  0.1790  1.0000   0.4449
# s(Temp0)         0        0.1143        0.0989  0.2647  0.3406   1.0000

# Temp0 concurved w Chl0
# removing Chl0

weekMod = gam(Pres ~ s(AEddyDist0,bs="ts",k=4)
              + s(CEddyDist0,bs="ts",k=4)
              # + s(Chl0,bs="ts",k=4)
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

# check convergence
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                 para s(AEddyDist0) s(CEddyDist0) s(Sal0) s(Temp0)
# para             1        0.0000        0.0000  0.0000   0.0000
# s(AEddyDist0)    0        1.0000        0.0058  0.0606   0.0234
# s(CEddyDist0)    0        0.0244        1.0000  0.1824   0.1164
# s(Sal0)          0        0.1933        0.1899  1.0000   0.4449
# s(Temp0)         0        0.1143        0.0989  0.3406   1.0000

# all <0.5

# Dredge all possible models nested within full model
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
#     fullSiteDayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
#                          + s(log(Chl0),bs="cs",k=5)
#                          + s(log(abs(FSLE0)),bs="cs",k=5)
#                          # + s(Sal0,bs="cs",k=5)
#                          + s(Sal700,bs="cs",k=5)
#                          + s(SSH0,bs="cs",k=5)
#                          + s(Temp0,bs="cs",k=4)
#                          # + s(Temp700,bs="cs",k=5)
#                          + s(sqrt(VelAsp0),bs="cs",k=5)
#                          + s(VelMag0,bs="cs",k=5),
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
#     png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'Daily.png',sep=""),width=600,height=600)
#     plot.gam(siteDayModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
#     while (dev.cur()>1) {dev.off()}
#     
#   }
#   
#   if (sum(weekData$Pres>0)>25){
#     #run full weekly model for this site
#     fullSiteWeekMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
#                           + s(log(Chl0),bs="cs",k=5)
#                           + s(log(abs(FSLE0)),bs="cs",k=5)
#                           # + s(Sal0,bs="cs",k=5)
#                           + s(Sal700,bs="cs",k=5)
#                           + s(SSH0,bs="cs",k=5)
#                           + s(Temp0,bs="cs",k=4)
#                           # + s(Temp700,bs="cs",k=5)
#                           + s(sqrt(VelAsp0),bs="cs",k=5)
#                           + s(VelMag0,bs="cs",k=5),
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
#     png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'Weekly.png',sep=""),width=600,height=600)
#     plot.gam(siteWeekModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
#     while (dev.cur()>1) {dev.off()}
#   }
# }
# 
# save(siteDayModList,pValDayList,siteDayModCompList,file=paste(outDir,'/',spec,'/','DailySiteSpecificModels.Rdata',sep=""))
# save(siteWeekModList,pValWeekList,siteWeekModCompList,file=paste(outDir,'/',spec,'/','WeeklySiteSpecificModels.Rdata',sep=""))
