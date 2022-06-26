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

## GAM approach ---------------------
# Regional model
spec = 'SpermWhale'
outDir = "J:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('J:/Chpt_3/ModelData/SpermWhale_masterDF.csv'))
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

# Remove zeros in FSLE data to prepare for later transformation
data$FSLE0[data$FSLE0==0] = NA
weeklyDF$FSLE0[weeklyDF$FSLE0==0] = NA

# Transform data to fix skew, get all predictors on a similar scale
weeklyDF$log_Chl0 = log10(weeklyDF$Chl0)
weeklyDF$log_abs_FSLE0 = log10(abs(weeklyDF$FSLE0))
weeklyDF$sqrt_CEddyDist0 = sqrt(weeklyDF$CEddyDist0)
weeklyDF$sqrt_AEddyDist0 = sqrt(weeklyDF$AEddyDist0)
weeklyDF$sqrt_VelAsp0 = sqrt(weeklyDF$VelAsp0)
weeklyDF$sqrt_VelAsp700 = sqrt(weeklyDF$VelAsp700)
weeklyDF$sqrt_EKE0 = sqrt(weeklyDF$EKE0)
weeklyDF$GSDist_div100 = weeklyDF$GSDist/100

# Remove incomplete observations (NAs in FSLE)
badRows = which(is.na(data),arr.ind=TRUE)[,1]
data = data[-badRows,]
badRows = unique(which(is.na(weeklyDF),arr.ind=TRUE)[,1])
weeklyDF = weeklyDF[-badRows,]

# re-round presence data
weeklyDF$Pres = round(weeklyDF$Pres)

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("sqrt_CEddyDist0",
                  "sqrt_EKE0",
                  "log_Chl0",
                  "Sal0",
                  "Sal700",
                  "SSH0",
                  "Temp0",
                  "Temp700")

# check residual autocorrelation of weekly data
sites = unique(weeklyDF$Site)
residAutocorr = matrix(ncol=1,nrow=length(sites))
rownames(residAutocorr) = sites
for (j in 1:length(sites)){
  
  siteInd = which(!is.na(str_match(weeklyDF$Site,sites[j])))
  
  if (sum(which(weeklyDF$Pres[siteInd]>0))>10){
    
    siteData = weeklyDF[siteInd,]
    
    BlockMod = glm(Pres~bs(sqrt_CEddyDist0) # include all terms in smoothVarList above!!
                   + bs(sqrt_EKE0)
                   + bs(log_Chl0)
                   + bs(Sal0)
                   + bs(Sal700)
                   + bs(SSH0)
                   + bs(Temp0)
                   + bs(Temp700),
                   data=siteData,family=poisson)
    
    acorr = acf(residuals(BlockMod), lag.max = 1000, main=paste(spec,"at",sites[j]))
    CI = ggfortify:::confint.acf(acorr)
    ACFidx = which(acorr[["acf"]] < CI, arr.ind=TRUE)
    residAutocorr[j,1] = ACFidx[1]
  }
}
residAutocorr
# HZ     5
# OC     3
# NC     3
# BC     3
# WC     4
# NFC    4
# HAT    3
# GS     2
# BP     2
# BS     2


# test for overdispersion
dispMod = glm(Pres~bs(sqrt_CEddyDist0) # include all terms in smoothVarList above!!
              + bs(sqrt_EKE0)
              + bs(log_Chl0)
              + bs(Sal0)
              + bs(Sal700)
              + bs(SSH0)
              + bs(Temp0)
              + bs(Temp700),
              data=weeklyDF,family=poisson)

dispersiontest(dispMod,alternative='two.sided')
# Dispersion test
# 
# data:  dispMod
# z = 15.262, p-value < 2.2e-16
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 145.2162 

# data are quite overdispersed, will use Tweedie family in models
modFam=tw

# determine whether covars should be included as linear or smooth terms
# seafloor aspect and water velocity direction are included as cyclic smooths

modOpts = c("linMod","threeKnots","fourKnots","fiveKnots")
AIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"Asp")){
    bs = "cc"
  } else { bs = "cs"}
  
  modelCall = paste('gam(Pres~weeklyDF$',smoothVarList[i],',data=weeklyDF,family=modFam)',sep="")
  linMod = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(weeklyDF$',smoothVarList[i],',bs="',bs,'",k=3),data=weeklyDF,family=modFam)',sep="")
  smoothMod1 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(weeklyDF$',smoothVarList[i],',bs="',bs,'",k=4),data=weeklyDF,family=modFam)',sep="")
  smoothMod2 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(weeklyDF$',smoothVarList[i],',bs="',bs,'",k=5),data=weeklyDF,family=modFam)',sep="")
  smoothMod3 = eval(parse(text=modelCall))
  
  AIC_votes[i,1:4] = c(AIC(linMod)[[1]],AIC(smoothMod1)[[1]],AIC(smoothMod2)[[1]],AIC(smoothMod3)[[1]])
  AIC_votes[i,5] = modOpts[which.min(AIC_votes[i,1:4])]
}

colnames(AIC_votes) = c(modOpts,"Best")
rownames(AIC_votes) = smoothVarList[]
AIC_votes
    
#                   linMod             threeKnots         fourKnots          fiveKnots          Best       
# sqrt_CEddyDist0 "18387.5392596658" "18374.4236498563" "18341.0614509327" "18343.6067063728" "fourKnots"
# sqrt_EKE0       "18469.1598271923" "18470.2594272735" "18470.2399536212" "18470.2007379425" "linMod"   
# log_Chl0        "18034.5133072321" "17989.2258390666" "17973.6788881557" "17956.1416829473" "fiveKnots"
# Sal0            "18175.0957277037" "18027.1706077071" "17918.9317834192" "17854.8549261408" "fiveKnots"
# Sal700          "18175.553955928"  "18092.7416136247" "18011.7382087124" "17986.8309561307" "fiveKnots"
# SSH0            "17755.5978432638" "17667.645068212"  "17668.5043198053" "17665.2386689044" "fiveKnots"
# Temp0           "18259.0806248987" "18260.8247905524" "18261.0078849554" "18247.7063269815" "fiveKnots"
# Temp700         "18036.4973610055" "18001.9671572893" "17946.1878704712" "17951.6835987591" "fourKnots"

# run full model
# fullMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
#               + s(log(Chl0),bs="cs",k=5)
#               + s(log(abs(FSLE0)),bs="cs",k=5)
#               + s(Sal0,bs="cs",k=4)
#               + s(Sal700,bs="cs",k=5)
#               + s(SSH0,bs="cs",k=5)
#               + s(Temp0,bs="cs",k=5)
#               + s(Temp700,bs="cs",k=5),
#               data=data,
#               family=poisson,
#               method="REML",
#               select=TRUE,
#               gamma=1.4,
#               na.action="na.fail")

# test whether to include EKE and EddyDist separately, or as an interaction term
g1 = gam(Pres~ s(sqrt_CEddyDist0,bs="cs",k=4)
         + sqrt_EKE0,
         data=weeklyDF,
         family=modFam,
         method="REML",
         select=TRUE,
         gamma=1.4,
         na.action="na.fail")

g2 = gam(Pres~ s(sqrt_CEddyDist0,sqrt_EKE0),
         data=weeklyDF,
         family=modFam,
         method="REML",
         select=TRUE,
         gamma=1.4,
         na.action="na.fail")

summary(g1)

# Family: Tweedie(p=1.785) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(sqrt_CEddyDist0, bs = "cs", k = 4) + sqrt_EKE0
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.175098   0.057056  90.702   <2e-16 ***
#   sqrt_EKE0   -0.003007   0.003169  -0.949    0.343    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df    F p-value    
# s(sqrt_CEddyDist0) 2.854      3 45.6  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0782   Deviance explained = 7.08%
# -REML = 6561.8  Scale est. = 4.0202    n = 1509

summary(g2)

# Family: Tweedie(p=1.787) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(sqrt_CEddyDist0, sqrt_EKE0)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.13585    0.02989   171.8   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df    F p-value    
# s(sqrt_CEddyDist0,sqrt_EKE0) 10.35     29 4.25  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0696   Deviance explained = 6.53%
# -REML = 6569.1  Scale est. = 4.0263    n = 1509

# slightly more deviance explained by including them separately

weekMod = gam(Pres ~ s(sqrt_CEddyDist0,bs="cs",k=4)
              + sqrt_EKE0
              + s(log_Chl0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=4),
              data=weeklyDF,
              family=modFam,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

# check convergence
# fullMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                       para s(sqrt_CEddyDist0) s(log_Chl0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700)
# para                  1             0.0000      0.0000  0.0000    0.0000  0.0000   0.0000     0.0000
# s(sqrt_CEddyDist0)    0             1.0000      0.0823  0.1817    0.1821  0.1010   0.0791     0.2945
# s(log_Chl0)           0             0.1707      1.0000  0.2348    0.2177  0.2120   0.2900     0.4315
# s(Sal0)               0             0.3502      0.2797  1.0000    0.6916  0.3354   0.2617     0.8075
# s(Sal700)             0             0.3492      0.3001  0.6968    1.0000  0.4132   0.2330     0.8535
# s(SSH0)               0             0.3117      0.3717  0.4210    0.4023  1.0000   0.3026     0.7603
# s(Temp0)              0             0.2078      0.3618  0.3551    0.3179  0.2631   1.0000     0.5035
# s(Temp700)            0             0.3309      0.2983  0.5394    0.5896  0.3448   0.3359     1.0000

# Sal0 concurved with Sal700,  Temp700
# Sal700 concurved with Sal0, Temp700
# Temp700 concurved with Sal0, Sal700
# removing Temp700, Sal0


# run reduced model
# dayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=4)
#              + s(log(Chl0),bs="cs",k=5)
#              + s(log(abs(FSLE0)),bs="cs",k=5)
#              # + s(Sal0,bs="cs",k=4)
#              + s(Sal700,bs="cs",k=5)
#              + s(SSH0,bs="cs",k=5)
#              + s(Temp0,bs="cs",k=5),
#              # + s(Temp700,bs="cs",k=5),
#              data=data,
#              family=poisson,
#              method="REML",
#              select=TRUE,
#              gamma=1.4,
#              na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt_CEddyDist0,bs="cs",k=4)
              + sqrt_EKE0
              + s(log_Chl0,bs="cs",k=5)
              # + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5),
              # + s(Temp700,bs="cs",k=4),
              data=weeklyDF,
              family=modFam,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

# check convergence
# dayMod$converged
# TRUE
weekMod$converged
# TRUE

conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                       para s(sqrt_CEddyDist0) s(log_Chl0) s(Sal700) s(SSH0) s(Temp0)
# para                  1             0.0000      0.0000    0.0000  0.0000   0.0000
# s(sqrt_CEddyDist0)    0             1.0000      0.0823    0.1821  0.1010   0.0791
# s(log_Chl0)           0             0.1707      1.0000    0.2177  0.2120   0.2900
# s(Sal700)             0             0.3492      0.3001    1.0000  0.4132   0.2330
# s(SSH0)               0             0.3117      0.3717    0.4023  1.0000   0.3026
# s(Temp0)              0             0.2078      0.3618    0.3179  0.2631   1.0000

# everything looks good! all values below 0.5

# dayModCompTable = dredge(dayMod,
#                          beta="none",
#                          evaluate=TRUE,
#                          trace=TRUE)

weekModCompTable = dredge(weekMod,
                          beta="none",
                          evaluate=TRUE,
                          trace=TRUE)

# run optimal models
# optDayMod = get.models(dayModCompTable,subset=1)
# optDayMod = optDayMod[[names(optDayMod)]]
# save(optDayMod,dayModCompTable,file=paste(outDir,'/',spec,'/','DailyRegionalModel.Rdata',sep=""))
optWeekMod = get.models(weekModCompTable,subset=1)
optWeekMod = optWeekMod[[names(optWeekMod)]]
save(optWeekMod,weekModCompTable,file=paste(outDir,'/',spec,'/','WeeklyRegionalModel.Rdata',sep=""))

sink(paste(outDir,'/',spec,'/','WeeklyRegionalModelSummary.txt',sep=""))
print(summary(optWeekMod))
sink()

# check p-values
# summary(optDayMod)

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(log(Chl0), 
#                                                   bs = "cs", k = 5) + s(Sal700, bs = "cs", k = 5) + 
#   s(sqrt(CEddyDist0), bs = "cs", k = 4) + s(SSH0, bs = "cs", 
#                                             k = 5) + s(Temp0, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 3.259379   0.002627    1241   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                       edf       Ref.df  Chi.sq  p-value    
#   s(log(abs(FSLE0)))  3.974      4      2412    <2e-16 ***
#   s(log(Chl0))        3.964      4      4217    <2e-16 ***
#   s(Sal700)           3.994      4      13736   <2e-16 ***
#   s(sqrt(CEddyDist0)) 2.997      3      3912    <2e-16 ***
#   s(SSH0)             3.996      4      22785   <2e-16 ***
#   s(Temp0)            3.991      4      15442   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.531   Deviance explained = 58.5%
# -REML =  86053  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: Tweedie(p=1.727) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log_Chl0, bs = "cs", k = 5) + s(Sal700, bs = "cs", k = 5) + 
#   s(sqrt_CEddyDist0, bs = "cs", k = 4) + s(SSH0, bs = "cs", 
#                                            k = 5) + s(Temp0, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.78312    0.02554   187.3   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df       F  p-value    
# s(log_Chl0)        0.5938      4   0.382  0.08788 .  
# s(Sal700)          0.9020      4   1.421  0.00706 ** 
#   s(sqrt_CEddyDist0) 2.6882      3   8.905 2.21e-06 ***
#   s(SSH0)            3.2873      4 126.632  < 2e-16 ***
#   s(Temp0)           3.6787      4  24.469  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.279   Deviance explained = 40.6%
# -REML = 6279.5  Scale est. = 3.4594    n = 1509

# plot
# png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesDaily.png',sep=""),width=600,height=600)
# plot.gam(optDayMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
# while (dev.cur()>1) {dev.off()}
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
    fullSiteDayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=4)
                         + s(log(Chl0),bs="cs",k=5)
                         + s(log(abs(FSLE0)),bs="cs",k=5)
                         # + s(Sal0,bs="cs",k=4)
                         + s(Sal700,bs="cs",k=5)
                         + s(SSH0,bs="cs",k=5)
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
    fullSiteWeekMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=4)
                          + s(log(Chl0),bs="cs",k=5)
                          + s(log(abs(FSLE0)),bs="cs",k=5)
                          # + s(Sal0,bs="cs",k=4)
                          + s(Sal700,bs="cs",k=5)
                          + s(SSH0,bs="cs",k=5)
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

