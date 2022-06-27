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

# library(tidyverse)
# library(splines2)
# library(geepack)
# source("getPvalues.R")

spec = 'Risso'
data = data.frame(read.csv('Risso_masterDF.csv'))
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
weeklyDF$sqrt_VelAsp400 = sqrt(weeklyDF$VelAsp400)
weeklyDF$sqrt_EKE0 = sqrt(weeklyDF$EKE0)
weeklyDF$GSDist_div100 = weeklyDF$GSDist/100

# Remove incomplete observations (NAs in FSLE)
badRows = which(is.na(data),arr.ind=TRUE)[,1]
data = data[-badRows,]
badRows = unique(which(is.na(weeklyDF),arr.ind=TRUE)[,1])
weeklyDF = weeklyDF[-badRows,]

# re-round presence data
weeklyDF$Pres = round(weeklyDF$Pres)


## GAM approach ---------------------
# Regional model
outDir = "J:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("sqrt_CEddyDist0",
                  "log_Chl0",
                  "sqrt_EKE0",
                  "log_abs_FSLE0",
                  "Sal0",
                  "Sal400",
                  "SSH0",
                  "Temp0",
                  "Temp400",
                  "VelMag400")

# check residual autocorrelation of weekly data
sites = unique(weeklyDF$Site)
residAutocorr = matrix(ncol=1,nrow=length(sites))
rownames(residAutocorr) = sites
for (j in 1:length(sites)){
  
  siteInd = which(!is.na(str_match(weeklyDF$Site,sites[j])))
  
  if (sum(which(weeklyDF$Pres[siteInd]>0))>10){
    
    siteData = weeklyDF[siteInd,]
    
    BlockMod = glm(Pres~ bs(sqrt_CEddyDist0) # include all terms in smoothVarList above!!
                   + bs(log_Chl0)
                   + bs(sqrt_EKE0)
                   + bs(log_abs_FSLE0)
                   + bs(Sal0)
                   + bs(Sal400)
                   + bs(SSH0)
                   + bs(Temp0)
                   + bs(Temp400)
                   + bs(VelMag400),
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
# NC     4
# BC     3
# WC     4
# NFC    7
# HAT    2
# GS     2
# BP     3
# BS     2

# test for overdispersion
dispMod = glm(Pres~ bs(sqrt_CEddyDist0) # include all terms in smoothVarList above!!
              + bs(log_Chl0)
              + bs(sqrt_EKE0)
              + bs(log_abs_FSLE0)
              + bs(Sal0)
              + bs(Sal400)
              + bs(SSH0)
              + bs(Temp0)
              + bs(Temp400)
              + bs(VelMag400),
              data=weeklyDF,family=poisson)

dispersiontest(dispMod,alternative='two.sided')
# Dispersion test
# 
# data:  dispMod
# z = 6.1553, p-value = 7.495e-10
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 103.5314  

# data are very overdispersed, will use Tweedie family in models

modFam=tw

# Test for how a term should be included in the model
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
# sqrt_CEddyDist0 "14358.3996101466" "14358.2835638967" "14342.1946394771" "14331.8187016489" "fiveKnots" 
# log_Chl0        "14318.4978623908" "14201.2941156678" "14186.6442252901" "14178.6019463055" "fiveKnots" 
# sqrt_EKE0       "14576.8784649556" "14574.581280182"  "14574.9512819411" "14575.0294498562" "threeKnots"
# log_abs_FSLE0   "14541.3849123199" "14400.3327352841" "14402.1422958548" "14370.979588982"  "fiveKnots" 
# Sal0            "14189.2031737403" "13644.7529892328" "13569.9752888658" "13573.4963366521" "fourKnots" 
# Sal400          "14170.2859121696" "13574.4350862965" "13501.6313090171" "13506.0806336159" "fourKnots" 
# SSH0            "13567.2205090723" "13568.7246205056" "13566.7945753462" "13566.226811799"  "fiveKnots" 
# Temp0           "14425.3413997134" "14220.5040598121" "14209.2518150626" "14212.5700477649" "fourKnots" 
# Temp400         "13973.3388184901" "13733.2212582477" "13725.4414969815" "13715.6388020514" "fiveKnots" 
# VelMag400       "14376.5276488689" "14376.2249115475" "14365.4266230264" "14353.7829734435" "fiveKnots" 

g1 = gam(Pres~ s(sqrt_CEddyDist0,bs="cs",k=5)
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
# Family: Tweedie(p=1.706) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(sqrt_CEddyDist0, bs = "cs", k = 5) + s(sqrt_EKE0, bs = "cs", 
#                                                   k = 3)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.00294    0.03692   108.4   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(sqrt_CEddyDist0) 3.278      4 74.31  < 2e-16 ***
#   s(sqrt_EKE0)       1.443      2 14.10 3.36e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.136   Deviance explained = 15.3%
# -REML = 5084.7  Scale est. = 6.5308    n = 1509

summary(g2)
# Family: Tweedie(p=1.706) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(sqrt_CEddyDist0, sqrt_EKE0)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.99777    0.03695   108.2   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(sqrt_CEddyDist0,sqrt_EKE0) 10.09     29 13.46  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.13   Deviance explained = 15.6%
# -REML = 5087.7  Scale est. = 6.521     n = 1509

# hardly any improvement in deviance explained by making an interaction, going to keep them separate
# run full model
# fullMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
#               + s(log(abs(FSLE0)),bs="cs",k=5)
#               + s(Sal400,bs="cs",k=5)
#               + s(SSH0,bs="cs",k=5)
#               + s(Temp0,bs="cs",k=5)
#               + s(VelMag400,bs="cc",k=5)
#               + s(sqrt(AEddyDist0),bs="cs",k=5)
#               + s(sqrt(CEddyDist0),bs="cs",k=5),
#               data=data,
#               family=poisson,
#               method="REML",
#               select=TRUE,
#               gamma=1.4,
#               na.action="na.fail")

weekMod = gam(Pres ~ + s(sqrt_CEddyDist0,bs="cs",k=5)
              + s(log_Chl0,bs="cs",k=5)
              + s(sqrt_EKE0,bs="cs",k=3)
              + s(log_abs_FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal400,bs="cs",k=4)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4)
              + s(Temp400,bs="cs",k=5)
              + s(VelMag400,bs="cc",k=5),
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

#                     para s(sqrt_CEddyDist0) s(log_Chl0) s(sqrt_EKE0) s(log_abs_FSLE0) s(Sal0) s(Sal400) s(SSH0) s(Temp0) s(Temp400) s(VelMag400)
# para                  1             0.0000      0.0000       0.0000           0.0000  0.0000    0.0000  0.0000   0.0000     0.0000       0.0000
# s(sqrt_CEddyDist0)    0             1.0000      0.0824       0.0506           0.0568  0.2333    0.2807  0.1073   0.1218     0.1597       0.1069
# s(log_Chl0)           0             0.0910      1.0000       0.0289           0.0758  0.2852    0.3192  0.2120   0.4581     0.2482       0.0574
# s(sqrt_EKE0)          0             0.0289      0.0127       1.0000           0.0204  0.0687    0.0780  0.0222   0.0126     0.0402       0.0358
# s(log_abs_FSLE0)      0             0.0728      0.0774       0.0462           1.0000  0.2195    0.2606  0.1792   0.1203     0.1782       0.1449
# s(Sal0)               0             0.1874      0.2743       0.0897           0.1512  1.0000    0.8966  0.2671   0.4006     0.4975       0.2013
# s(Sal400)             0             0.1942      0.2819       0.0913           0.1586  0.8737    1.0000  0.2696   0.3521     0.4984       0.2161
# s(SSH0)               0             0.1804      0.3717       0.0789           0.1620  0.5051    0.5658  1.0000   0.4799     0.4606       0.3286
# s(Temp0)              0             0.1114      0.3568       0.0209           0.0858  0.3514    0.3730  0.2578   1.0000     0.4618       0.1354
# s(Temp400)            0             0.1778      0.3042       0.0701           0.1688  0.6160    0.6930  0.3257   0.5376     1.0000       0.2112
# s(VelMag400)          0             0.1028      0.0596       0.0964           0.0861  0.2721    0.3160  0.2940   0.1439     0.1706       1.0000
# > 

# Sal0 concurved w Sal400, SSH, Temp400,
# Sal400 concurved w Sal0, SSH, Temp400
# Temp0 concurved w Temp400

# taking out Sal0, Temp400, VelMag400

# dayMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
#               + s(log(abs(FSLE0)),bs="cs",k=5)
#               + s(Sal400,bs="cs",k=5)
#               # + s(SSH0,bs="cs",k=5)
#               + s(Temp0,bs="cs",k=5)
#               + s(VelMag400,bs="cc",k=5)
#               + s(sqrt(AEddyDist0),bs="cs",k=5)
#               + s(sqrt(CEddyDist0),bs="cs",k=5),
#               data=data,
#               family=poisson,
#               method="REML",
#               select=TRUE,
#               gamma=1.4,
#               na.action="na.fail")

weekMod = gam(Pres ~ + s(sqrt_CEddyDist0,bs="cs",k=5)
              + s(log_Chl0,bs="cs",k=5)
              + s(sqrt_EKE0,bs="cs",k=3)
              + s(log_abs_FSLE0,bs="cs",k=5)
              # + s(Sal0,bs="cs",k=4)
              + s(Sal400,bs="cs",k=4)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4),
              # + s(Temp400,bs="cs",k=5)
              # + s(VelMag400,bs="cc",k=5),
              data=weeklyDF,
              family=modFam,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                      para s(sqrt_CEddyDist0) s(log_Chl0) s(sqrt_EKE0) s(log_abs_FSLE0) s(Sal400) s(SSH0) s(Temp0)
# para                  1             0.0000      0.0000       0.0000           0.0000    0.0000  0.0000   0.0000
# s(sqrt_CEddyDist0)    0             1.0000      0.0824       0.0506           0.0568    0.2807  0.1073   0.1218
# s(log_Chl0)           0             0.0910      1.0000       0.0289           0.0758    0.3192  0.2120   0.4581
# s(sqrt_EKE0)          0             0.0289      0.0127       1.0000           0.0204    0.0780  0.0222   0.0126
# s(log_abs_FSLE0)      0             0.0728      0.0774       0.0462           1.0000    0.2606  0.1792   0.1203
# s(Sal400)             0             0.1942      0.2819       0.0913           0.1586    1.0000  0.2696   0.3521
# s(SSH0)               0             0.1804      0.3717       0.0789           0.1620    0.5658  1.0000   0.4799
# s(Temp0)              0             0.1114      0.3568       0.0209           0.0858    0.3730  0.2578   1.0000

# taking out SSH

weekMod = gam(Pres ~ + s(sqrt_CEddyDist0,bs="cs",k=5)
              + s(log_Chl0,bs="cs",k=5)
              + s(sqrt_EKE0,bs="cs",k=3)
              + s(log_abs_FSLE0,bs="cs",k=5)
              # + s(Sal0,bs="cs",k=4)
              + s(Sal400,bs="cs",k=4)
              # + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4),
              # + s(Temp400,bs="cs",k=5)
              # + s(VelMag400,bs="cc",k=5),
              data=weeklyDF,
              family=modFam,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                     para s(sqrt_CEddyDist0) s(log_Chl0) s(sqrt_EKE0) s(log_abs_FSLE0) s(Sal400) s(Temp0)
# para                  1             0.0000      0.0000       0.0000           0.0000    0.0000   0.0000
# s(sqrt_CEddyDist0)    0             1.0000      0.0824       0.0506           0.0568    0.2807   0.1218
# s(log_Chl0)           0             0.0910      1.0000       0.0289           0.0758    0.3192   0.4581
# s(sqrt_EKE0)          0             0.0289      0.0127       1.0000           0.0204    0.0780   0.0126
# s(log_abs_FSLE0)      0             0.0728      0.0774       0.0462           1.0000    0.2606   0.1203
# s(Sal400)             0             0.1942      0.2819       0.0913           0.1586    1.0000   0.3521
# s(Temp0)              0             0.1114      0.3568       0.0209           0.0858    0.3730   1.0000

# no convurvity >0.5

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
#                                                   k = 5) + s(Sal400, bs = "cs", k = 5) + s(sqrt(AEddyDist0), 
#                                                                                            bs = "cs", k = 5) + s(sqrt(CEddyDist0), bs = "cs", k = 5) + 
#   s(Temp0, bs = "cs", k = 5) + s(VelMag400, bs = "cc", k = 5) + 
#   1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.38674    0.00865   160.3   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df   Chi.sq p-value    
# s(log(abs(FSLE0)))  3.976      4  2986.93  <2e-16 ***
#   s(log(Chl0))        3.985      4  2800.66  <2e-16 ***
#   s(Sal400)           3.998      4 14959.60  <2e-16 ***
#   s(sqrt(AEddyDist0)) 3.932      4   454.93  <2e-16 ***
#   s(sqrt(CEddyDist0)) 3.827      4  1371.22  <2e-16 ***
#   s(Temp0)            3.987      4  1080.52  <2e-16 ***
#   s(VelMag400)        2.419      3    91.89  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.208   Deviance explained = 37.8%
# -REML =  64062  Scale est. = 1         n = 10375

summary(optWeekMod)

#Family: Tweedie(p=1.641) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log_abs_FSLE0, bs = "cs", k = 5) + s(log_Chl0, bs = "cs", 
#                                                 k = 5) + s(Sal400, bs = "cs", k = 4) + s(sqrt_CEddyDist0, 
#                                                                                          bs = "cs", k = 5) + s(sqrt_EKE0, bs = "cs", k = 3) + s(Temp0, 
#                                                                                                                                                 bs = "cs", k = 4) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.37554    0.03447   97.94   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df       F  p-value    
# s(log_abs_FSLE0)   3.280      4  14.959 6.34e-14 ***
#   s(log_Chl0)        3.408      4   8.021 1.71e-07 ***
#   s(Sal400)          2.970      3 189.359  < 2e-16 ***
#   s(sqrt_CEddyDist0) 1.004      4   2.063  0.00216 ** 
#   s(sqrt_EKE0)       1.408      2   8.853 1.46e-05 ***
#   s(Temp0)           2.341      3   3.797  0.00289 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.295   Deviance explained = 49.2%
# -REML = 4780.5  Scale est. = 5.0838    n = 1509

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
                         + s(Sal400,bs="cs",k=5)
                         + s(Temp0,bs="cs",k=5)
                         + s(VelMag400,bs="cc",k=5)
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
    fullSiteWeekMod = gam(Pres ~ s(log(Chl0),bs="cs",k=5)
                          + s(log(abs(FSLE0)),bs="cs",k=5)
                          + s(Sal400,bs="cs",k=5)
                          + s(Temp0,bs="cs",k=5)
                          + s(VelMag400,bs="cc",k=5)
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





# GEEGLM approach ------------------------------------
lagID = 43
numClust = length(data$Pres)/(lagID-1)
if (numClust<length(data$Pres)){
  clustID = rep(1:ceiling(numClust),each=lagID)
  clustID = clustID[1:length(data$Pres)]
} else {
  clustID = 1:length(data$Pres)
}
data$GroupID = clustID

# Round presence to get Poisson dist
data$Pres = round(data$Pres)

# Test for how a term should be included in the model
startTime = Sys.time()
smoothVarList = c("Temp0",
                  "Sal0",
                  "Sal200",
                  "SSH0",
                  "FSLE0",
                  "GSLat",
                  "Slope",
                  "Aspect")

modOpts = c("linMod","threeKnots","fourKnots")
QIC_votes = matrix(nrow=length(smoothVarList),ncol=4)

for (i in 1:(length(smoothVarList)-1)){
  
  modelCall = paste('geeglm(Pres~data$',smoothVarList[i],',data=data,family=poisson,id=GroupID,corstr="ar1")',sep="")
  linMod = eval(parse(text=modelCall))
  
  modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.5)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],'))),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  smoothMod1 = eval(parse(text=modelCall))
  
  modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.333,0.666)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],'))),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  smoothMod2 = eval(parse(text=modelCall))
  
  # modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.275,0.5,0.725)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],'))),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  # smoothMod3 = eval(parse(text=modelCall))
  
  QIC_votes[i,1:3] = c(QIC(linMod)[[1]],QIC(smoothMod1)[[1]],QIC(smoothMod2)[[1]])
  QIC_votes[i,4] = modOpts[which.min(QIC_votes[i,1:3])]
}

endTime = Sys.time()
endTime-startTime
colnames(QIC_votes) = c(modOpts,"Best")
rownames(QIC_votes) = smoothVarList[]
QIC_votes

#           linMod              threeKnots          fourKnots           Best        
# Temp0  "-2458133.44127198" "-2497958.76610649" "-2494447.98398555" "threeKnots"
# Sal0   "-2387145.82259727" "-2397529.44160543" "-2397770.58164437" "fourKnots" 
# Sal200 "-2385452.46919174" "-2396045.24119207" "-2397783.78049046" "fourKnots" 
# SSH0   "-2509932.72715967" "-2510614.47948344" "-2512094.69825831" "fourKnots" 
# FSLE0  "-2385328.82959341" "-2386058.84470483" "-2386081.53007312" "fourKnots" 
# GSLat  "-2386298.16752796" "-2387378.29830821" "-2387372.85097507" "threeKnots"
# Slope  "-2404596.46021069" "-2530688.29829264" "-2532163.73313361" "fourKnots" 
# Aspect NA                  NA                  NA                  NA          

# Make smooth terms, run full model and check collinearity
smoothVarList = c("Temp0",
                  "Sal0",
                  "Sal200",
                  "SSH0",
                  "FSLE0",
                  "GSLat",
                  "Slope",
                  "Aspect")
knotList = list(c(0.5),
                c(0.333,0.666),
                c(0.333,0.666),
                c(0.333,0.666),
                c(0.333,0.666),
                c(0.5),
                c(0.333,0.666),
                c(0.333,0.666))
linVarList = list()
smoothNameList = character()


for (i in 1:length(smoothVarList)){
  
  if (str_detect(smoothVarList[i],"Asp")){
    eval(parse(text=paste('S_',smoothVarList[i],'= mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=unlist(knotList[i])),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')),periodic=TRUE)',sep="")))
  } else {
    eval(parse(text=paste('S_',smoothVarList[i],'= mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=unlist(knotList[i])),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')))',sep="")))
  }
  
  smoothNameList = c(smoothNameList,paste('S_',smoothVarList[i],sep=""))
}
thisForm = formula(paste('Pres~',paste(c(smoothNameList,linVarList),collapse="+"),sep=""))

fullMod = geeglm(thisForm,
                 family=poisson,
                 data=data,
                 id=GroupID,
                 corstr="unstructured")

VIFvals = vif(fullMod)
VIFvals = cbind(VIFvals,(VIFvals[,3])^2)
colnames(VIFvals)[4] = "LOOK AT ME"
VIFvals

#               GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# S_Temp0   8.585876  4        1.308347   1.711773
# S_Sal0   31.063539  5        1.410019   1.988155
# S_Sal200 34.562502  5        1.425150   2.031052
# S_SSH0   49.168005  5        1.476278   2.179398
# S_FSLE0   4.956267  5        1.173587   1.377308
# S_GSLat   2.698404  4        1.132109   1.281672
# S_Slope  88.670559  5        1.565950   2.452200
# S_Aspect 11.408530  2        1.837839   3.377652

# no collinearity, no need to remove any terms

# check convergence
fullMod$geese$error
# 1 model didn't converge

# run w. unstructured correlation
fullMod = geeglm(thisForm,
                 family=poisson,
                 data=data,
                 id=GroupID,
                 corstr="unstructured")
fullMod$geese$error
# 0 model converged

#               GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# S_Temp0   5.896675  4        1.248320   1.558303
# S_Sal0   14.987336  5        1.310909   1.718482
# S_Sal200 16.682475  5        1.325031   1.755707
# S_SSH0   30.147910  5        1.405807   1.976293
# S_FSLE0   3.262515  5        1.125525   1.266807
# S_GSLat   2.855038  4        1.140122   1.299879
# S_Slope  78.268675  5        1.546532   2.391760
# S_Aspect 12.849682  2        1.893316   3.584645

# check term significance
PV = getPvalues(fullMod)

#Remove non-significant terms, re-run model, check p-values
PV$'p-value'[PV$'p-value'=="<0.0001"] = 0.0001
badVars = PV$Variable[as.numeric(PV$'p-value')>=0.05]
redMod<-eval(parse(text=paste("update(fullMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
if (redMod$geese$error==1){
  print("Model did not converge")
} else {PVred = getPvalues(redMod)
PVred$'p-value'[PVred$'p-value'=="<0.0001"] = 0.0001
PVred}

# Get p-values
PV = getPvalues(reducedMod)

# Plot terms from regional model
source("plotSmooths.R")
source("plotLinears.R")
terms = names(reducedMod$model)[2:length(names(reducedMod$model))]
k=3

for (i in 1:length(terms)){
  if (str_detect(terms[i],"S_")){ # plot smooth terms
    term = str_remove(terms[i],"S_")
    coefInd = which(str_detect(names(reducedMod$coefficients),term))
    if (str_detect(term,"Asp")){periodic=TRUE} else {periodic=FALSE}
    print(plotSmooths(reducedMod,term,coefInd,k,periodic,site=NA,title=NULL))
  } else { # plot linear terms
    term=terms[i]
    coefInd = which(str_detect(names(reducedMod$coefficients),term))
    print(plotLinears(reducedMod,term,coefInd,site=NA,title=NULL))
  }
}
