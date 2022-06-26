library(tidyverse)
library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)
library(forecast)
library(nlme)
library(itsadug)
library(AER)

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
weeklyDF$GSDist_div100 = weeklyDF$GSDist/100

# Remove incomplete observations (NAs in FSLE)
badRows = which(is.na(data),arr.ind=TRUE)[,1]
data = data[-badRows,]
badRows = unique(which(is.na(weeklyDF),arr.ind=TRUE)[,1])
weeklyDF = weeklyDF[-badRows,]

# re-round presence data
weeklyDF$Pres = round(weeklyDF$Pres)

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c ("sqrt_CEddyDist0",
                   "EKE0",
                   "log_abs_FSLE0",
                   "Sal0",
                   "Sal700",
                   "SSH0",
                   "Temp0",
                   "Temp700",
                   "sqrt_VelAsp0",
                   "VelMag0")

# check residual autocorrelation of weekly data
sites = unique(weeklyDF$Site)
residAutocorr = matrix(ncol=1,nrow=length(sites))
rownames(residAutocorr) = sites
for (j in 1:length(sites)){
  
  siteInd = which(!is.na(str_match(weeklyDF$Site,sites[j])))
  
  if (sum(which(weeklyDF$Pres[siteInd]>0))>10){
    
    siteData = weeklyDF[siteInd,]
    
    BlockMod = glm(Pres~bs(sqrt_CEddyDist0) # include all terms in smoothVarList above!!
                   + bs(EKE0)
                   + bs(log_abs_FSLE0)
                   + bs(Sal0)
                   + bs(Sal700)
                   + bs(SSH0)
                   + bs(Temp0) 
                   + bs(Temp700)
                   + bs(sqrt_VelAsp0)
                   + bs(VelMag0),
                   data=siteData,family=poisson)
    
    acorr = acf(residuals(BlockMod), lag.max = 1000, main=paste(spec,"at",sites[j]))
    CI = ggfortify:::confint.acf(acorr)
    ACFidx = which(acorr[["acf"]] < CI, arr.ind=TRUE)
    residAutocorr[j,1] = ACFidx[1]
  }
}
residAutocorr

# HZ    NA
# OC    NA
# NC    NA
# BC    NA
# WC    NA
# NFC   NA
# HAT    2
# GS     2
# BP     2
# BS     4

# test for overdispersion
dispMod = glm(Pres~bs(sqrt_CEddyDist0) # include all terms in smoothVarList above!!
              + bs(EKE0)
              + bs(log_abs_FSLE0)
              + bs(Sal0)
              + bs(Sal700)
              + bs(SSH0)
              + bs(Temp0) 
              + bs(Temp700)
              + bs(sqrt_VelAsp0)
              + bs(VelMag0),
              data=weeklyDF,family=poisson)

dispersiontest(dispMod,alternative='two.sided')
# Dispersion test
# 
# data:  dispMod
# z = 7.4955, p-value = 6.605e-14
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 5.873926  
# data are somewhat overdispersed, will use Tweedie family in models
modFam=tw

# Test for how a term should be included in the model
modOpts = c("linMod","threeKnots","fourKnots","fiveKnots")
AIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"Asp")){
    bs = "cc"
  } else { bs = "cs"}
  
  modelCall = paste('gam(Pres~weeklyDF$',smoothVarList[i],',data=weeklyDF,family=modFam,method="REML",select=TRUE)',sep="")
  linMod = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(weeklyDF$',smoothVarList[i],',bs="',bs,'",k=3),data=weeklyDF,family=modFam,method="REML",select=TRUE)',sep="")
  smoothMod1 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(weeklyDF$',smoothVarList[i],',bs="',bs,'",k=4),data=weeklyDF,family=modFam,method="REML",select=TRUE)',sep="")
  smoothMod2 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(weeklyDF$',smoothVarList[i],',bs="',bs,'",k=5),data=weeklyDF,family=modFam,method="REML",select=TRUE)',sep="")
  smoothMod3 = eval(parse(text=modelCall))
  
  AIC_votes[i,1:4] = c(AIC(linMod)[[1]],AIC(smoothMod1)[[1]],AIC(smoothMod2)[[1]],AIC(smoothMod3)[[1]])
  AIC_votes[i,5] = modOpts[which.min(AIC_votes[i,1:4])]
}

colnames(AIC_votes) = c(modOpts,"Best")
rownames(AIC_votes) = smoothVarList[]
AIC_votes

#                   linMod             threeKnots         fourKnots          fiveKnots          Best        
# sqrt_CEddyDist0 "3219.90565013409" "3219.12816702401" "3220.03331829053" "3220.16671593889" "threeKnots"
# EKE0            "3199.94911972609" "3199.85470391814" "3201.4535400976"  "3199.84856165046" "fiveKnots" 
# log_abs_FSLE0   "3245.42599611214" "3236.50565357833" "3212.74081960895" "3210.48532122779" "fiveKnots" 
# Sal0            "2673.64016800021" "2670.25928342262" "2670.29494820878" "2671.46470888367" "threeKnots"
# Sal700          "2520.68207020349" "2444.07858113014" "2444.11273748937" "2443.70237307537" "fiveKnots" 
# SSH0            "2670.52653164919" "2586.10389773372" "2586.28275418063" "2588.68647898231" "threeKnots"
# Temp0           "3128.91818565504" "3022.22649582107" "2993.24172873932" "3001.98376449282" "fourKnots" 
# Temp700         "2985.22043420561" "2781.16366681897" "2785.06864153591" "2772.03920547426" "fiveKnots" 
# sqrt_VelAsp0    "3245.83316696644" "3251.40037117783" "3251.40037117783" "3250.9154076727"  "linMod"    
# VelMag0         "3213.42905771662" "3198.92444196861" "3197.61866504164" "3198.52963159969" "fourKnots"

# run full model
# fullMod = gam(Pres ~ s(sqrt(EKE0),bs="cs",k=5)
#               + s(Sal0,bs="cs",k=3)
#               + s(Sal700,bs="cs",k=5)
#               + s(Temp0,bs="cs",k=5)
#               + s(sqrt(VelAsp700),bs="cc",k=5)
#               + s(VelMag700,bs="cs",k=5)
#               + s(sqrt(AEddyDist0),bs="cs",k=5)
#               + s(sqrt(CEddyDist0),bs="cs",k=5),
#               data=data,
#               family=poisson,
#               gamma=1.4,
#               na.action="na.fail",
#               method="REML",
#               select=TRUE)

weekMod = gam(Pres ~ s(sqrt_CEddyDist0,bs="cs",k=3)
              + s(EKE0,bs="cs",k=5)
              + s(log_abs_FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=3)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=3)
              + s(Temp0,bs="cs",k=4)
              + s(Temp700,bs="cs",k=5)
              + s(sqrt_VelAsp0,bs="cc",k=4) # including as smooth so it can be cyclic
              + s(VelMag0,bs="cs",k=4),
              data=weeklyDF,
              family=modFam,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              select=TRUE)

# check convergence
# fullMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                     para s(sqrt_CEddyDist0) s(EKE0) s(log_abs_FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(sqrt_VelAsp0) s(VelMag0)
# para                  1             0.0000  0.0000           0.0000  0.0000    0.0000  0.0000   0.0000     0.0000          0.0000     0.0000
# s(sqrt_CEddyDist0)    0             1.0000  0.1098           0.0400  0.2896    0.1683  0.2268   0.1114     0.1547          0.1233     0.1306
# s(EKE0)               0             0.1230  1.0000           0.0858  0.2468    0.2130  0.2534   0.1310     0.1663          0.3845     0.9362
# s(log_abs_FSLE0)      0             0.0901  0.1857           1.0000  0.2661    0.1782  0.2585   0.1203     0.1835          0.2398     0.1895
# s(Sal0)               0             0.2610  0.2077           0.1450  1.0000    0.5643  0.6764   0.3869     0.5511          0.3494     0.2884
# s(Sal700)             0             0.2718  0.4396           0.1787  0.8823    1.0000  0.7638   0.3853     0.6364          0.4498     0.4862
# s(SSH0)               0             0.2093  0.3325           0.1457  0.7104    0.3393  1.0000   0.4687     0.4159          0.3675     0.3281
# s(Temp0)              0             0.1605  0.1744           0.0858  0.4405    0.3104  0.5783   1.0000     0.3825          0.2410     0.2321
# s(Temp700)            0             0.2631  0.2645           0.1672  0.8523    0.6399  0.7551   0.4374     1.0000          0.3983     0.3568
# s(sqrt_VelAsp0)       0             0.1254  0.4177           0.1337  0.3917    0.2575  0.3907   0.1858     0.2298          1.0000     0.5099
# s(VelMag0)            0             0.1274  0.8301           0.0876  0.2602    0.2192  0.2655   0.1370     0.1739          0.3934     1.0000

# EKE0 problematic with VelMag0
# Sal0 problematic w Sal700, SSH, Temp700
# Sal700 problematic w Sal0, Temp700
# SSH concurved w Sal0, Sal700, Temp700
# Temp700 concurved w Sal0, Sal700
# VelMag700 problematic w EKE0, Sal700, VelAsp700

# All equally significant, arbitrarily taking out Sal0, VelMag700, Temp700, SSH

# run reduced model
# dayMod = gam(Pres ~ 
#                # s(sqrt(EKE0),bs="cs",k=5)
#              # + s(Sal0,bs="cs",k=3)
#              + s(Sal700,bs="cs",k=5)
#              + s(Temp0,bs="cs",k=5)
#              + s(sqrt(VelAsp700),bs="cc",k=5)
#              + s(VelMag700,bs="cs",k=5)
#              + s(sqrt(AEddyDist0),bs="cs",k=5)
#              + s(sqrt(CEddyDist0),bs="cs",k=5),
#               data=data,
#               family=poisson,
#               gamma=1.4,
#               na.action="na.fail",
#               method="REML",
#               select=TRUE)

weekMod = gam(Pres ~ s(sqrt_CEddyDist0,bs="cs",k=3)
              + s(EKE0,bs="cs",k=5)
              + s(log_abs_FSLE0,bs="cs",k=5)
              # + s(Sal0,bs="cs",k=3)
              + s(Sal700,bs="cs",k=5)
              # + s(SSH0,bs="cs",k=3)
              + s(Temp0,bs="cs",k=4)
              # + s(Temp700,bs="cs",k=5)
              + s(sqrt_VelAsp0,bs="cc",k=4),
              # + s(VelMag0,bs="cs",k=4),
             data=weeklyDF,
             family=modFam,
             gamma=1.4,
             na.action="na.fail",
             method="REML",
             select=TRUE)

# check convergence
# dayMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                     para s(sqrt_CEddyDist0) s(EKE0) s(log_abs_FSLE0) s(Sal700) s(Temp0) s(sqrt_VelAsp0)
# para                  1             0.0000  0.0000           0.0000    0.0000   0.0000          0.0000
# s(sqrt_CEddyDist0)    0             1.0000  0.1098           0.0400    0.1683   0.1114          0.1233
# s(EKE0)               0             0.1230  1.0000           0.0858    0.2130   0.1310          0.3845
# s(log_abs_FSLE0)      0             0.0901  0.1857           1.0000    0.1782   0.1203          0.2398
# s(Sal700)             0             0.2718  0.4396           0.1787    1.0000   0.3853          0.4498
# s(Temp0)              0             0.1605  0.1744           0.0858    0.3104   1.0000          0.2410
# s(sqrt_VelAsp0)       0             0.1254  0.4177           0.1337    0.2575   0.1858          1.0000

# only slight concurvity remains

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
# 
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

# Family: Tweedie(p=1.475) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(EKE0, bs = "cs", k = 5) + s(log_abs_FSLE0, bs = "cs", 
#                                        k = 5) + s(Sal700, bs = "cs", k = 5) + s(sqrt_VelAsp0, bs = "cc", 
#                                                                                 k = 4) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   -50.16      23.35  -2.148   0.0319 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
# s(EKE0)          1.731      4  8.572  < 2e-16 ***
#   s(log_abs_FSLE0) 2.204      4  4.114 0.000144 ***
#   s(Sal700)        1.921      4 31.858  < 2e-16 ***
#   s(sqrt_VelAsp0)  1.753      2 11.129 3.72e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.524   Deviance explained = 83.9%
# -REML = 853.87  Scale est. = 8.1659    n = 1509

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
