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
spec = 'Cuvier'
outDir = "E:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('Cuvier_masterDF.csv'))
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
smoothVarList = c ("sqrt_CEddyDist0",
                   "sqrt_EKE0",
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
    
    BlockMod = glm(Pres~bs(sqrt_CEddyDist0)
                   + bs(sqrt_EKE0)
                   + bs(log_abs_FSLE0) # include all terms in smoothVarList above!!
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
# HZ     3
# OC     3
# NC     2
# BC     3
# WC     4
# NFC    2
# HAT    3
# GS     3
# BP    NA
# BS     2


# test for overdispersion
dispMod = glm(Pres~bs(sqrt_CEddyDist0)
              + bs(sqrt_EKE0)
              + bs(log_abs_FSLE0) # include all terms in smoothVarList above!!
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
# z = 1.4116, p-value = 0.1581
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 1351932 

# data are incredibly overdispersed, will use Tweedie family in models
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
# sqrt_CEddyDist0 "9754.72423406743" "9756.26467145516" "9756.47148128688" "9756.45235198655" "linMod"    
# sqrt_EKE0       "9860.83422089882" "9851.7002937517"  "9852.08746501125" "9852.1636702777"  "threeKnots"
# log_abs_FSLE0   "9696.46004906687" "9604.72947978713" "9606.25225797437" "9594.59946081962" "fiveKnots" 
# Sal0            "9743.50130199551" "9737.05998812076" "9583.78162202874" "9454.40483636207" "fiveKnots" 
# Sal700          "9765.95731096871" "9735.56729801862" "9537.43367247707" "9395.99999615024" "fiveKnots" 
# SSH0            "9879.64029705624" "9264.56785077703" "9245.09687783305" "8851.71788310579" "fiveKnots" 
# Temp0           "9811.71759192573" "9813.01054681113" "9798.76837192849" "9801.13951303969" "fourKnots" 
# Temp700         "9757.16519239658" "9754.29972609094" "9603.67884182524" "9606.21502137614" "fourKnots" 
# sqrt_VelAsp0    "9497.35938508613" "9510.04362168551" "9510.04362168551" "9485.38447162852" "fiveKnots" 
# VelMag0         "9237.35017636557" "9239.14185260165" "9226.16617605049" "9217.7568665327"  "fiveKnots" 

g1 = gam(Pres~ sqrt_CEddyDist0
         + s(sqrt_EKE0,bs="cs",k=3),
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
# Family: Tweedie(p=1.722) 
# Link function: log 
# 
# Formula:
#   Pres ~ sqrt_CEddyDist0 + s(sqrt_EKE0, bs = "cs", k = 3)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      5.55041    0.22607   24.55   <2e-16 ***
#   sqrt_CEddyDist0 -0.20827    0.01887  -11.04   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(sqrt_EKE0) 1.613      2 10.13 6.86e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0514   Deviance explained = 9.02%
# -REML = 3477.8  Scale est. = 10.262    n = 1509

summary(g2)

# run full model
# fullMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
#               + s(log(abs(FSLE0)),bs="cs",k=5)
#               + s(Sal0,bs="cs",k=5)
#               + s(Sal700,bs="cs",k=5)
#               + s(SSH0,bs="cs",k=5)
#               + s(Temp0,bs="cs",k=5)
#               + s(Temp700,bs="cs",k=4),
#               data=data,
#               family=poisson,
#               method="REML",
#               select=TRUE,
#               gamma=1.4,
#               na.action="na.fail")

weekMod = gam(Pres ~ sqrt_CEddyDist0
              + s(sqrt_EKE0,bs="cs",k=3)
              + s(log_abs_FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4)
              + s(Temp700,bs="cs",k=4)
              + s(sqrt_VelAsp0,bs="cc",k=5)
              + s(VelMag0,bs="cs",k=5),
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

#                   para s(sqrt_EKE0) s(log_abs_FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(sqrt_VelAsp0) s(VelMag0)
# para                1       0.0000           0.0000  0.0000    0.0000  0.0000   0.0000     0.0000          0.0000     0.0000
# s(sqrt_EKE0)        0       1.0000           0.0204  0.0499    0.0443  0.0222   0.0126     0.0682          0.0220     0.0423
# s(log_abs_FSLE0)    0       0.0462           1.0000  0.1767    0.1782  0.1792   0.1203     0.3143          0.1217     0.1418
# s(Sal0)             0       0.0899           0.1720  1.0000    0.6916  0.3354   0.4311     0.8075          0.2144     0.3392
# s(Sal700)           0       0.0827           0.1787  0.6968    1.0000  0.4132   0.3853     0.8535          0.2252     0.3869
# s(SSH0)             0       0.0789           0.1620  0.4210    0.4023  1.0000   0.4799     0.7603          0.2178     0.3791
# s(Temp0)            0       0.0209           0.0858  0.3375    0.3104  0.2578   1.0000     0.5001          0.1179     0.1776
# s(Temp700)          0       0.0778           0.1657  0.5394    0.5896  0.3448   0.4229     1.0000          0.2004     0.2537
# s(sqrt_VelAsp0)     0       0.0448           0.1349  0.2709    0.2636  0.2367   0.1890     0.4340          1.0000     0.3952
# s(VelMag0)          0       0.1472           0.0878  0.2257    0.2204  0.3237   0.1374     0.3172          0.2079     1.0000

# Sal0 concurved with Sal700 and Temp700
# Sal700 concurved with Sal0 and Temp700
# Temp700 concurved with Sal700, Sal0, SSH, Temp0  

# removing Sal0, Temp700, VelMag0


# run reduced model
# dayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
#              + s(log(abs(FSLE0)),bs="cs",k=5)
#              # + s(Sal0,bs="cs",k=5)
#              + s(Sal700,bs="cs",k=5)
#              + s(SSH0,bs="cs",k=5)
#              + s(Temp0,bs="cs",k=5),
#              # + s(Temp700,bs="cs",k=4),
#              data=data,
#              family=poisson,
#              method="REML",
#              select=TRUE,
#              gamma=1.4,
#              na.action="na.fail")

weekMod = gam(Pres ~ sqrt_CEddyDist0
              + s(sqrt_EKE0,bs="cs",k=3)
              + s(log_abs_FSLE0,bs="cs",k=5)
              # + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4)
              # + s(Temp700,bs="cs",k=4)
              + s(sqrt_VelAsp0,bs="cc",k=5),
              # + s(VelMag0,bs="cs",k=5),
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

#                   para s(sqrt_EKE0) s(log_abs_FSLE0) s(Sal700) s(SSH0) s(Temp0) s(sqrt_VelAsp0)
# para                1       0.0000           0.0000    0.0000  0.0000   0.0000          0.0000
# s(sqrt_EKE0)        0       1.0000           0.0204    0.0443  0.0222   0.0126          0.0220
# s(log_abs_FSLE0)    0       0.0462           1.0000    0.1782  0.1792   0.1203          0.1217
# s(Sal700)           0       0.0827           0.1787    1.0000  0.4132   0.3853          0.2252
# s(SSH0)             0       0.0789           0.1620    0.4023  1.0000   0.4799          0.2178
# s(Temp0)            0       0.0209           0.0858    0.3104  0.2578   1.0000          0.1179
# s(sqrt_VelAsp0)     0       0.0448           0.1349    0.2636  0.2367   0.1890          1.0000

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
# save(optWeekMod,weekModCompTable,file=paste(outDir,'/',spec,'/','WeeklyRegionalModel.Rdata',sep=""))
save(optWeekMod,weekModCompTable,file=paste(spec,'_','WeeklyRegionalModel_Updated.Rdata',sep=""))

# sink(paste(outDir,'/',spec,'/','WeeklyRegionalModelSummary.txt',sep=""))
sink(paste(spec,'_','WeeklyRegionalModelSummary_Updated.txt',sep=""))
print(summary(optWeekMod))
sink()

# check p-values
# summary(optDayMod)

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

# Family: Tweedie(p=1.557) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(Sal700, bs = "cs", k = 5) + s(sqrt_VelAsp0, bs = "cc", 
#                                          k = 5) + s(SSH0, bs = "cs", k = 5) + s(Temp0, bs = "cs", 
#                                                                                 k = 4) + sqrt_CEddyDist0 + 1
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      2.80937    0.19208  14.626  < 2e-16 ***
#   sqrt_CEddyDist0 -0.08100    0.01561  -5.187 2.43e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F p-value    
# s(Sal700)       3.100      4  34.00  <2e-16 ***
#   s(sqrt_VelAsp0) 2.248      3  10.72  <2e-16 ***
#   s(SSH0)         3.942      4 277.79  <2e-16 ***
#   s(Temp0)        2.643      3  13.81  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.462   Deviance explained = 61.7%
# -REML = 3080.2  Scale est. = 6.1933    n = 1509

# plot
# png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesDaily.png',sep=""),width=600,height=600)
# plot.gam(optDayMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
# while (dev.cur()>1) {dev.off()}
# png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesWeekly.png',sep=""),width=600,height=600)
png(filename=paste(spec,'_allSitesWeekly_Updated.png',sep=""),width=600,height=600)
plot.gam(optWeekMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
while (dev.cur()>1) {dev.off()}

# Retrain model with 2/3 of data, then validate by predicting remaining 1/3
trainInd = sample(1:nrow(weeklyDF),floor(nrow(weeklyDF)*.66))
testInd = setdiff(1:nrow(weeklyDF),trainInd)

valMod = update(optWeekMod,data=weeklyDF[trainInd,])
plot.gam(valMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
save(valMod,weekModCompTable,file=paste(spec,'_','ValidationModel.Rdata',sep=""))

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

