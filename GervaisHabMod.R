library(tidyverse)
library(mgcv.helper)
library(splines)
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
spec = 'Gervais'
outDir = "E:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('Gervais_masterDF.csv'))
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
smoothVarList = c ("log_Chl0",
                   "sqrt_EKE0",
                   "log_abs_FSLE0",
                   "Sal0",
                   "Sal700",
                   "SSH0",
                   "Temp0",
                   "Temp700",
                   "sqrt_VelAsp0")

# check residual autocorrelation of weekly data
sites = unique(weeklyDF$Site)
residAutocorr = matrix(ncol=1,nrow=length(sites))
rownames(residAutocorr) = sites
for (j in 1:length(sites)){
  
  siteInd = which(!is.na(str_match(weeklyDF$Site,sites[j])))
  
  if (sum(which(weeklyDF$Pres[siteInd]>0))>10){
    
    siteData = weeklyDF[siteInd,]
    
    BlockMod = glm(Pres~bs(log_Chl0) # include all terms in smoothVarList above!!
                   + bs(sqrt_EKE0)
                   + bs(log_abs_FSLE0)
                   + bs(Sal0)
                   + bs(Sal700)
                   + bs(SSH0)
                   + bs(Temp0) 
                   + bs(Temp700)
                   + bs(sqrt_VelAsp0),
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
# NFC    2
# HAT    4
# GS     3
# BP     4
# BS     3


# test for overdispersion
dispMod = glm(Pres~bs(log_Chl0) # include all terms in smoothVarList above!!
              + bs(sqrt_EKE0)
              + bs(log_abs_FSLE0)
              + bs(Sal0)
              + bs(Sal700)
              + bs(SSH0)
              + bs(Temp0) 
              + bs(Temp700)
              + bs(sqrt_VelAsp0),
              data=weeklyDF,family=poisson)

dispersiontest(dispMod,alternative='two.sided')
# Dispersion test
# 
# data:  dispMod
# z = 1.2013, p-value = 0.2296
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 1.943917e+42  
# model didn't converge; will use Tweedie family in models anyway
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

#               linMod             threeKnots         fourKnots          fiveKnots          Best        
# log_Chl0      "7728.54856893595" "7476.49604097709" "7372.65571565143" "7366.20345530275" "fiveKnots" 
# sqrt_EKE0     "8086.93495104378" "8084.72494189301" "8085.15797667064" "8085.37608857843" "threeKnots"
# log_abs_FSLE0 "8016.246608491"   "8014.24113854719" "7919.15540948349" "7920.96360273379" "fourKnots" 
# Sal0          "6900.52000295644" "6903.42568994993" "6658.4059071418"  "6588.95349089539" "fiveKnots" 
# Sal700        "7007.27575457946" "6865.11799815084" "6665.95020371157" "6559.04840828499" "fiveKnots" 
# SSH0          "7346.29003683508" "6542.42798591582" "6518.89624221537" "6511.90997912821" "fiveKnots" 
# Temp0         "7457.56564439985" "7222.27867912965" "7217.32974886274" "7208.08642394956" "fiveKnots" 
# Temp700       "6744.4809590541"  "6602.30758392683" "6569.81133158467" "6575.20798858821" "fourKnots" 
# sqrt_VelAsp0  "7605.99835051968" "7555.61291532273" "7555.61291532273" "7541.18767942448" "fiveKnots" 

# run full model
# fullMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
#               + s(log(abs(FSLE0)),bs="cs",k=5)
#               + s(Sal0,bs="cs",k=4)
#               + s(Sal700,bs="cs",k=5)
#               + s(SSH0,bs="cs",k=5)
#               + s(Temp700,bs="cs",k=5),
#               data=data,
#               family=poisson,
#               gamma=1.4,
#               na.action="na.fail",
#               method="REML",
#               select=TRUE)

weekMod = gam(Pres ~ s(log_Chl0,bs="cs",k=5)
              + s(sqrt_EKE0,bs="cs",k=3)
              + s(log_abs_FSLE0,bs="cs",k=4)
              + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=4)
              + s(sqrt_VelAsp0,bs="cc",k=4),
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

#                    para s(log_Chl0) s(sqrt_EKE0) s(log_abs_FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(sqrt_VelAsp0)
# para                1      0.0000       0.0000           0.0000  0.0000    0.0000  0.0000   0.0000     0.0000          0.0000
# s(log_Chl0)         0      1.0000       0.0289           0.1426  0.2348    0.2177  0.2120   0.2900     0.4315          0.1652
# s(sqrt_EKE0)        0      0.0127       1.0000           0.0380  0.0499    0.0443  0.0222   0.0097     0.0682          0.0321
# s(log_abs_FSLE0)    0      0.0724       0.0461           1.0000  0.1736    0.1756  0.1754   0.0631     0.3099          0.2331
# s(Sal0)             0      0.2797       0.0899           0.3021  1.0000    0.6916  0.3354   0.2617     0.8075          0.4266
# s(Sal700)           0      0.3001       0.0827           0.3083  0.6968    1.0000  0.4132   0.2330     0.8535          0.4498
# s(SSH0)             0      0.3717       0.0789           0.2870  0.4210    0.4023  1.0000   0.3026     0.7603          0.4116
# s(Temp0)            0      0.3618       0.0209           0.1420  0.3551    0.3179  0.2631   1.0000     0.5035          0.2447
# s(Temp700)          0      0.2983       0.0778           0.3077  0.5394    0.5896  0.3448   0.3359     1.0000          0.3951
# s(sqrt_VelAsp0)     0      0.0999       0.0414           0.2401  0.2640    0.2575  0.2045   0.1186     0.4219          1.0000

# Sal0 problematic w Sal 700, Temp700
# Sal700 problematic w Temp700, Sal0
# Taking out Temp700, Sal0

# run reduced model
# dayMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
#              + s(log(abs(FSLE0)),bs="cs",k=5)
#              # + s(Sal0,bs="cs",k=4)
#              + s(Sal700,bs="cs",k=5)
#              + s(SSH0,bs="cs",k=5),
#              # + s(Temp700,bs="cs",k=5),
#              data=data,
#              family=poisson,
#              gamma=1.4,
#              na.action="na.fail",
#              method="REML",
#              select=TRUE)

weekMod = gam(Pres ~ s(log_Chl0,bs="cs",k=5)
              + s(sqrt_EKE0,bs="cs",k=3)
              + s(log_abs_FSLE0,bs="cs",k=4)
              # + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              # + s(Temp700,bs="cs",k=4)
              + s(sqrt_VelAsp0,bs="cc",k=4),
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

#                     para s(log_Chl0) s(sqrt_EKE0) s(log_abs_FSLE0) s(Sal700) s(SSH0) s(Temp0) s(sqrt_VelAsp0)
# para                1      0.0000       0.0000           0.0000    0.0000  0.0000   0.0000          0.0000
# s(log_Chl0)         0      1.0000       0.0289           0.1426    0.2177  0.2120   0.2900          0.1652
# s(sqrt_EKE0)        0      0.0127       1.0000           0.0380    0.0443  0.0222   0.0097          0.0321
# s(log_abs_FSLE0)    0      0.0724       0.0461           1.0000    0.1756  0.1754   0.0631          0.2331
# s(Sal700)           0      0.3001       0.0827           0.3083    1.0000  0.4132   0.2330          0.4498
# s(SSH0)             0      0.3717       0.0789           0.2870    0.4023  1.0000   0.3026          0.4116
# s(Temp0)            0      0.3618       0.0209           0.1420    0.3179  0.2631   1.0000          0.2447
# s(sqrt_VelAsp0)     0      0.0999       0.0414           0.2401    0.2575  0.2045   0.1186          1.0000

# all <0.5

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
# save(optWeekMod,weekModCompTable,file=paste(outDir,'/',spec,'/','WeeklyRegionalModel.Rdata',sep=""))
save(optWeekMod,weekModCompTable,file=paste(spec,'_','WeeklyRegionalModel.Rdata',sep=""))

# sink(paste(outDir,'/',spec,'/','WeeklyRegionalModelSummary.txt',sep=""))
sink(paste(spec,'_','WeeklyRegionalModelSummary.txt',sep=""))
print(summary(optWeekMod))
sink()

# check p-values
# summary(optDayMod)

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(log(Chl0), bs = "cs", 
#                                                   k = 5) + s(Sal700, bs = "cs", k = 5) + s(SSH0, bs = "cs", 
#                                                                                            k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.3783     0.5177  -8.457   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq p-value    
# s(log(abs(FSLE0))) 3.336      4   62.15  <2e-16 ***
#   s(log(Chl0))       3.982      4 1976.90  <2e-16 ***
#   s(Sal700)          3.624      4 6175.78  <2e-16 ***
#   s(SSH0)            3.040      4 4113.23  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.418   Deviance explained = 66.9%
# -REML =  22060  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: Tweedie(p=1.493) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log_Chl0, bs = "cs", k = 5) + s(Sal700, bs = "cs", k = 5) + 
#   s(sqrt_VelAsp0, bs = "cc", k = 4) + s(SSH0, bs = "cs", k = 5) + 
#   s(Temp0, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -3.1365     0.3737  -8.393   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
# s(log_Chl0)     3.2285      4 12.345 1.48e-11 ***
#   s(Sal700)       3.6791      4 42.887  < 2e-16 ***
#   s(sqrt_VelAsp0) 1.5833      2  5.152  0.00149 ** 
#   s(SSH0)         2.7814      4 37.214  < 2e-16 ***
#   s(Temp0)        0.6522      4  0.513  0.06545 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =    0.6   Deviance explained = 86.9%
# -REML = 2169.5  Scale est. = 5.3331    n = 1509

# plot
# png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesDaily.png',sep=""),width=600,height=600)
# plot.gam(optDayMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
# while (dev.cur()>1) {dev.off()}
# png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesWeekly.png',sep=""),width=600,height=600)
png(filename=paste(spec,'_allSitesWeekly.png',sep=""),width=600,height=600)
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
    fullSiteDayMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
                         + s(log(abs(FSLE0)),bs="cs",k=5)
                         + s(Sal700,bs="cs",k=5)
                         + s(SSH0,bs="cs",k=5),
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
    fullSiteWeekMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
                          + s(log(abs(FSLE0)),bs="cs",k=5)
                          + s(Sal700,bs="cs",k=5)
                          + s(SSH0,bs="cs",k=5),
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

