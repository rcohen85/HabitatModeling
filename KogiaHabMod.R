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
spec = 'Kogia'
outDir = "J:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('Kogia_masterDF.csv'))
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

# don't use observations from northern sites, could be harbor porpoise
southSites = weeklyDF$Site%in%c("HAT","GS","BP","BS")
weeklyDF = weeklyDF[southSites,]

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c ("sqrt_AEddyDist0",
                   "sqrt_CEddyDist0",
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
    
    BlockMod = glm(Pres~bs(sqrt_AEddyDist0) # include all terms in smoothVarList above!!
                   + bs(sqrt_CEddyDist0)
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
# HAT    2
# GS     3
# BP     3
# BS     2


# test for overdispersion
dispMod = glm(Pres~bs(sqrt_AEddyDist0) # include all terms in smoothVarList above!!
              + bs(sqrt_CEddyDist0)
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
# z = 1.7466, p-value = 0.0807
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 10.98851    

# data are overdispersed, will use Tweedie family in models
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
# sqrt_AEddyDist0 "3105.68863419029" "3078.87197980631" "3071.4329284973"  "3071.98014454789" "fourKnots" 
# sqrt_CEddyDist0 "3177.35375082314" "3177.3566473534"  "3177.32635574486" "3177.22313050125" "fiveKnots" 
# sqrt_EKE0       "3179.12849355677" "3177.21618937477" "3177.22006167014" "3177.21598981078" "fiveKnots" 
# log_Chl0        "3102.76644585721" "3046.10060550227" "3001.99775363892" "3003.50584077327" "fourKnots" 
# Sal0            "3033.7336703498"  "3017.06651610142" "3010.43445968513" "3011.31173002615" "fourKnots" 
# Sal700          "3044.54790211274" "3014.19375155709" "3015.83745771975" "3015.67152391201" "threeKnots"
# SSH0            "3049.24161953275" "2943.37214077998" "2914.35624265369" "2919.81973907783" "fourKnots" 
# Temp0           "3171.98881771714" "3165.92604996109" "3144.66070654918" "3141.18337926811" "fiveKnots" 
# Temp700         "3178.98009449532" "3106.93459940847" "3097.82167557185" "3096.81687346604" "fiveKnots" 

g1 = gam(Pres~ s(sqrt_AEddyDist0,bs="cs",k=4)
         + s(sqrt_CEddyDist0,bs="cs",k=5)
         + s(sqrt_EKE0,bs="cs",k=5),
         data=weeklyDF,
         family=modFam,
         method="REML",
         select=TRUE,
         gamma=1.4,
         na.action="na.fail")

g2 = gam(Pres~ s(sqrt_CEddyDist0,sqrt_EKE0)
         + s(sqrt_AEddyDist0,sqrt_EKE0),
         data=weeklyDF,
         family=modFam,
         method="REML",
         select=TRUE,
         gamma=1.4,
         na.action="na.fail")

summary(g1)
# Family: Tweedie(p=1.372) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(sqrt_AEddyDist0, bs = "cs", k = 4) + s(sqrt_CEddyDist0, 
#                                                   bs = "cs", k = 5) + s(sqrt_EKE0, bs = "cs", k = 5)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.3156     0.0489   26.91   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F p-value    
# s(sqrt_AEddyDist0) 2.634598      3 33.845  <2e-16 ***
#   s(sqrt_CEddyDist0) 0.503488      4  0.344  0.0948 .  
# s(sqrt_EKE0)       0.001125      4  0.000  1.0000    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0784   Deviance explained = 13.6%
# -REML = 1094.5  Scale est. = 3.1424    n = 615

summary(g2)
# Family: Tweedie(p=1.388) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(sqrt_CEddyDist0, sqrt_EKE0) + s(sqrt_AEddyDist0, sqrt_EKE0)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.36027    0.04866   27.95   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(sqrt_CEddyDist0,sqrt_EKE0) 0.5592     29 0.038   0.138    
# s(sqrt_AEddyDist0,sqrt_EKE0) 0.9833     28 2.879  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0235   Deviance explained = 9.09%
# -REML = 1103.9  Scale est. = 3.2436    n = 615

# more deviance explained by keeping them separate

# run full model
# fullMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
#               + s(log(Chl0),bs="cs",k=5)
#               + s(log(abs(FSLE0)),bs="cs",k=5)
#               + s(Sal0,bs="cs",k=5)
#               + s(Sal700,bs="cs",k=5)
#               + s(SSH0,bs="cs",k=5)
#               + s(Temp0,bs="cs",k=4)
#               + s(Temp700,bs="cs",k=5)
#               + s(sqrt(VelAsp0),bs="cs",k=5)
#               + s(VelMag0,bs="cs",k=5),
#               data=data,
#               family=poisson,
#               method="REML",
#               select=TRUE,
#               gamma=1.4,
#               na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt_AEddyDist0,bs="cs",k=4)
              + s(sqrt_CEddyDist0,bs="cs",k=5)
              + s(sqrt_EKE0,bs="cs",k=5)
              + s(log_Chl0,bs="cs",k=4)
              + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=3)
              + s(SSH0,bs="cs",k=4)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=5),
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

#                    para s(sqrt_AEddyDist0) s(sqrt_CEddyDist0) s(sqrt_EKE0) s(log_Chl0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700)
# para                  1             0.0000             0.0000       0.0000      0.0000  0.0000    0.0000  0.0000   0.0000     0.0000
# s(sqrt_AEddyDist0)    0             1.0000             0.0166       0.0019      0.1275  0.1406    0.1506  0.2353   0.0037     0.0439
# s(sqrt_CEddyDist0)    0             0.0238             1.0000       0.0091      0.0168  0.0301    0.0256  0.0767   0.0175     0.0179
# s(sqrt_EKE0)          0             0.0048             0.0085       1.0000      0.0118  0.0065    0.0060  0.0076   0.0127     0.0224
# s(log_Chl0)           0             0.1592             0.0264       0.0080      1.0000  0.1978    0.2736  0.3833   0.1204     0.1435
# s(Sal0)               0             0.1250             0.0160       0.0120      0.0480  1.0000    0.5501  0.2228   0.2062     0.0932
# s(Sal700)             0             0.1186             0.0028       0.0212      0.1765  0.3964    1.0000  0.3419   0.0209     0.0704
# s(SSH0)               0             0.2595             0.0363       0.0220      0.3980  0.3070    0.3992  1.0000   0.0737     0.2624
# s(Temp0)              0             0.0201             0.0161       0.0196      0.2106  0.3103    0.0642  0.0533   1.0000     0.3261
# s(Temp700)            0             0.0577             0.0242       0.0272      0.1406  0.1427    0.1530  0.1213   0.3033     1.0000

# Sal700 concurved w Sal0

# removing Sal700


# run reduced model
# dayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
#              + s(log(Chl0),bs="cs",k=5)
#              + s(log(abs(FSLE0)),bs="cs",k=5)
#              # + s(Sal0,bs="cs",k=5)
#              + s(Sal700,bs="cs",k=5)
#              + s(SSH0,bs="cs",k=5)
#              + s(Temp0,bs="cs",k=4)
#              # + s(Temp700,bs="cs",k=5)
#              + s(sqrt(VelAsp0),bs="cs",k=5)
#              + s(VelMag0,bs="cs",k=5),
#              data=data,
#              family=poisson,
#              method="REML",
#              select=TRUE,
#              gamma=1.4,
#              na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt_AEddyDist0,bs="cs",k=4)
              + s(sqrt_CEddyDist0,bs="cs",k=5)
              + s(sqrt_EKE0,bs="cs",k=5)
              + s(log_Chl0,bs="cs",k=4)
              + s(Sal0,bs="cs",k=4)
              # + s(Sal700,bs="cs",k=3)
              + s(SSH0,bs="cs",k=4)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=5),
              data=weeklyDF,
              family=modFam,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

# check convergence
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                    para s(sqrt_AEddyDist0) s(sqrt_CEddyDist0) s(sqrt_EKE0) s(log_Chl0) s(Sal0) s(SSH0) s(Temp0) s(Temp700)
# para                  1             0.0000             0.0000       0.0000      0.0000  0.0000  0.0000   0.0000     0.0000
# s(sqrt_AEddyDist0)    0             1.0000             0.0166       0.0019      0.1275  0.1406  0.2353   0.0037     0.0439
# s(sqrt_CEddyDist0)    0             0.0238             1.0000       0.0091      0.0168  0.0301  0.0767   0.0175     0.0179
# s(sqrt_EKE0)          0             0.0048             0.0085       1.0000      0.0118  0.0065  0.0076   0.0127     0.0224
# s(log_Chl0)           0             0.1592             0.0264       0.0080      1.0000  0.1978  0.3833   0.1204     0.1435
# s(Sal0)               0             0.1250             0.0160       0.0120      0.0480  1.0000  0.2228   0.2062     0.0932
# s(SSH0)               0             0.2595             0.0363       0.0220      0.3980  0.3070  1.0000   0.0737     0.2624
# s(Temp0)              0             0.0201             0.0161       0.0196      0.2106  0.3103  0.0533   1.0000     0.3261
# s(Temp700)            0             0.0577             0.0242       0.0272      0.1406  0.1427  0.1213   0.3033     1.0000

# no concurvity >0.5

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
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(log(Chl0), 
#                                                   bs = "cs", k = 5) + s(Sal700, bs = "cs", k = 5) + 
#   s(sqrt(CEddyDist0), bs = "cs", k = 5) + s(sqrt(VelAsp0), 
#                                             bs = "cs", k = 5) + s(SSH0, bs = "cs", k = 5) + 
#   s(Temp0, bs = "cs", k = 4) + s(VelMag0, bs = "cs", 
#                                  k = 5) + 1
# 
# Parametric coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.12169    0.03812  -55.66   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                        edf Ref.df  Chi.sq  p-value    
#   s(log(abs(FSLE0)))  3.1820      4  27.931 1.33e-06 ***
#   s(log(Chl0))        3.8406      4 123.049  < 2e-16 ***
#   s(Sal700)           0.9284      4   4.282   0.0149 *  
#   s(sqrt(CEddyDist0)) 2.9235      4  55.528  < 2e-16 ***
#   s(sqrt(VelAsp0))    2.2423      4   8.063   0.0146 *  
#   s(SSH0)             3.2541      4 221.357  < 2e-16 ***
#   s(Temp0)            2.7327      3 156.350  < 2e-16 ***
#   s(VelMag0)          3.1196      4  44.798  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.0794   Deviance explained = 31.1%
# -REML = 4623.8  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: Tweedie(p=1.316) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log_Chl0, bs = "cs", k = 4) + s(Sal0, bs = "cs", k = 4) + 
#   s(sqrt_AEddyDist0, bs = "cs", k = 4) + s(sqrt_CEddyDist0, 
#                                            bs = "cs", k = 5) + s(SSH0, bs = "cs", k = 4) + s(Temp700, 
#                                                                                              bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.90240    0.06055    14.9   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
# s(log_Chl0)        0.6959      3  0.777   0.0525 .  
# s(Sal0)            2.3378      3  6.422 2.05e-05 ***
#   s(sqrt_AEddyDist0) 2.0583      3  6.370 2.40e-05 ***
#   s(sqrt_CEddyDist0) 0.7079      4  0.673   0.0458 *  
#   s(SSH0)            2.8838      3 25.037  < 2e-16 ***
#   s(Temp700)         1.6765      4  6.323 2.24e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.185   Deviance explained = 40.4%
# -REML =   1018  Scale est. = 2.5006    n = 615

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
    fullSiteDayMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=5)
                         + s(log(Chl0),bs="cs",k=5)
                         + s(log(abs(FSLE0)),bs="cs",k=5)
                         # + s(Sal0,bs="cs",k=5)
                         + s(Sal700,bs="cs",k=5)
                         + s(SSH0,bs="cs",k=5)
                         + s(Temp0,bs="cs",k=4)
                         # + s(Temp700,bs="cs",k=5)
                         + s(sqrt(VelAsp0),bs="cs",k=5)
                         + s(VelMag0,bs="cs",k=5),
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
                          + s(log(Chl0),bs="cs",k=5)
                          + s(log(abs(FSLE0)),bs="cs",k=5)
                          # + s(Sal0,bs="cs",k=5)
                          + s(Sal700,bs="cs",k=5)
                          + s(SSH0,bs="cs",k=5)
                          + s(Temp0,bs="cs",k=4)
                          # + s(Temp700,bs="cs",k=5)
                          + s(sqrt(VelAsp0),bs="cs",k=5)
                          + s(VelMag0,bs="cs",k=5),
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
