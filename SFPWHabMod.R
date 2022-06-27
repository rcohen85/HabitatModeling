library(tidyverse)
library(mgcv.helper)
library(splines2)
library(splines)
library(mgcv)
library(MuMIn)
library(gratia)
library(forecast)
library(nlme)
library(itsadug)
library(AER)

## GAM approach ---------------------
# Regional model
spec = 'SFPW'
outDir = "J:/Chpt_3/GAM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('UD26_masterDF.csv'))
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
smoothVarList = c ("sqrt_AEddyDist0",
                   "sqrt_CEddyDist0",
                   "sqrt_EKE0",
                   "log_Chl0",
                   "Sal0",
                   "Sal400",
                   "SSH0",
                   "Temp0",
                   "Temp400")

# check residual autocorrelation of weekly data
sites = unique(weeklyDF$Site)
residAutocorr = matrix(ncol=1,nrow=length(sites))
rownames(residAutocorr) = sites
for (j in 1:length(sites)){
  
  siteInd = which(!is.na(str_match(weeklyDF$Site,sites[j])))
  
  if (sum(which(weeklyDF$Pres[siteInd]>0))>10){
    
    siteData = weeklyDF[siteInd,]
    
    BlockMod = glm(Pres~bs(sqrt_CEddyDist0)
                   + bs(sqrt_AEddyDist0) # include all terms in smoothVarList above!!
                   + bs(sqrt_EKE0)
                   + bs(log_Chl0)
                   + bs(Sal0)
                   + bs(Sal400)
                   + bs(SSH0)
                   + bs(Temp0)
                   + bs(Temp400),
                   data=siteData,family=poisson)
    
    acorr = acf(residuals(BlockMod), lag.max = 1000, main=paste(spec,"at",sites[j]))
    CI = ggfortify:::confint.acf(acorr)
    ACFidx = which(acorr[["acf"]] < CI, arr.ind=TRUE)
    residAutocorr[j,1] = ACFidx[1]
  }
}
residAutocorr

# HZ     2
# OC     3
# NC     2
# BC     4
# WC     2
# NFC    6
# HAT    3
# GS     2
# BP     2
# BS     2

# test for overdispersion
dispMod = glm(Pres~bs(sqrt_CEddyDist0)
              + bs(sqrt_AEddyDist0) # include all terms in smoothVarList above!!
              + bs(sqrt_EKE0)
              + bs(log_Chl0)
              + bs(Sal0)
              + bs(Sal400)
              + bs(SSH0)
              + bs(Temp0)
              + bs(Temp400),
              data=weeklyDF,family=poisson)

dispersiontest(dispMod,alternative='two.sided')
# Dispersion test
# 
# data:  dispMod
# z = 10.308, p-value < 2.2e-16
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 41.04422 

# data are pretty overdispersed, will use Tweedie family in models
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
# sqrt_AEddyDist0 "12722.582201656"  "12721.572143237"  "12658.634262172"  "12663.862428509"  "fourKnots" 
# sqrt_CEddyDist0 "13008.4264609154" "13007.4044317204" "13007.5199293844" "13007.4667491667" "threeKnots"
# sqrt_EKE0       "13044.4372678205" "13044.7819753503" "13044.7905733883" "13044.7817405871" "linMod"    
# log_Chl0        "12775.0389369442" "12605.2701656248" "12607.4852259379" "12608.394711936"  "threeKnots"
# Sal0            "13043.2455317858" "13033.7068730274" "12924.0386507564" "12880.4978528141" "fiveKnots" 
# Sal400          "13047.7862860471" "12989.5473330095" "12889.0896065478" "12771.2785862178" "fiveKnots" 
# SSH0            "12809.777793419"  "12215.0056026009" "12216.7745682014" "12209.2116757598" "fiveKnots" 
# Temp0           "13043.1663572304" "12913.3341690512" "12895.9720229951" "12888.9533402921" "fiveKnots" 
# Temp400         "13047.4883458692" "12927.7799674183" "12866.6402070648" "12866.2759576736" "fiveKnots"  

# test whether to include EKE and EddyDist separately, or as an interaction term
g1 = gam(Pres~ s(sqrt_CEddyDist0,bs="cs",k=3)
         + s(sqrt_AEddyDist0,bs="cs",k=4)
         + sqrt_EKE0,
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
# Family: Tweedie(p=1.687) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(sqrt_CEddyDist0, bs = "cs", k = 3) + s(sqrt_AEddyDist0, 
#                                                   bs = "cs", k = 4) + sqrt_EKE0
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 3.179897   0.070516  45.095  < 2e-16 ***
#   sqrt_EKE0   0.010052   0.003803   2.644  0.00829 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(sqrt_CEddyDist0) 1.435      2  21.2 1.21e-11 ***
#   s(sqrt_AEddyDist0) 2.921      3 153.2  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.207   Deviance explained = 20.7%
# -REML = 4485.3  Scale est. = 5.8312    n = 1509

summary(g2)
# Family: Tweedie(p=1.69) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(sqrt_CEddyDist0, sqrt_EKE0) + s(sqrt_AEddyDist0, sqrt_EKE0)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.34582    0.03761   88.97   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
# s(sqrt_CEddyDist0,sqrt_EKE0)  1.389     29  1.586 3.03e-12 ***
#   s(sqrt_AEddyDist0,sqrt_EKE0) 13.873     28 15.689  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.19   Deviance explained = 20.2%
# -REML = 4498.6  Scale est. = 5.8609    n = 1509

# not much difference, keeping them separate

# run full model
# fullMod = gam(Pres ~ s(sqrt(CEddyDist0),bs="cs",k=4)
#               + s(log(Chl0),bs="cs",k=5)
#               + s(log(abs(FSLE0)),bs="cs",k=5)
#               # + s(Sal0,bs="cs",k=4)
#               + s(Sal700,bs="cs",k=5)
#               + s(SSH0,bs="cs",k=5)
#               + s(Temp0,bs="cs",k=5),
#               # + s(Temp700,bs="cs",k=5),
#               data=data,
#               family=poisson,
#               method="REML",
#               select=TRUE,
#               gamma=1.4,
#               na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt_AEddyDist0,bs="cs",k=4)
              + s(sqrt_CEddyDist0,bs="cs",k=3)
              + sqrt_EKE0
              + s(log_Chl0,bs="cs",k=3)
              + s(Sal0,bs="cs",k=5)
              + s(Sal400,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp400,bs="cs",k=5),
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

#                       para s(sqrt_AEddyDist0) s(sqrt_CEddyDist0) s(log_Chl0) s(Sal0) s(Sal400) s(SSH0) s(Temp0) s(Temp400)
# para                  1             0.0000             0.0000      0.0000  0.0000    0.0000  0.0000   0.0000     0.0000
# s(sqrt_AEddyDist0)    0             1.0000             0.0084      0.0852  0.0547    0.0770  0.0626   0.0284     0.0721
# s(sqrt_CEddyDist0)    0             0.0068             1.0000      0.1107  0.1703    0.1782  0.0942   0.0719     0.1499
# s(log_Chl0)           0             0.1199             0.1162      1.0000  0.2124    0.2043  0.1825   0.2430     0.2048
# s(Sal0)               0             0.1563             0.2699      0.3667  1.0000    0.8392  0.3354   0.2617     0.5212
# s(Sal400)             0             0.2298             0.2755      0.3821  0.8332    1.0000  0.3830   0.2443     0.5214
# s(SSH0)               0             0.1850             0.2373      0.5282  0.4210    0.4327  1.0000   0.3026     0.4606
# s(Temp0)              0             0.0871             0.1630      0.5032  0.3551    0.3263  0.2631   1.0000     0.4850
# s(Temp400)            0             0.1581             0.2551      0.4366  0.5396    0.5597  0.3257   0.4776     1.0000

# Chl0 concurved w SSH, Temp0
# Sal0 concurved with Sal400 and Temp400
# Sal400 concurved with Sal0 and Temp400
# Temp400 concurved with Sal0, Sal400

# remove Temp400, Sal0, SSH


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

weekMod = gam(Pres ~ s(sqrt_AEddyDist0,bs="cs",k=4)
              + s(sqrt_CEddyDist0,bs="cs",k=3)
              + sqrt_EKE0
              + s(log_Chl0,bs="cs",k=3)
              # + s(Sal0,bs="cs",k=5)
              + s(Sal400,bs="cs",k=5)
              # + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5),
              # + s(Temp400,bs="cs",k=5),
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

#                     para s(sqrt_AEddyDist0) s(sqrt_CEddyDist0) s(log_Chl0) s(Sal400) s(Temp0)
# para                  1             0.0000             0.0000      0.0000    0.0000   0.0000
# s(sqrt_AEddyDist0)    0             1.0000             0.0084      0.0852    0.0770   0.0284
# s(sqrt_CEddyDist0)    0             0.0068             1.0000      0.1107    0.1782   0.0719
# s(log_Chl0)           0             0.1199             0.1162      1.0000    0.2043   0.2430
# s(Sal400)             0             0.2298             0.2755      0.3821    1.0000   0.2443
# s(Temp0)              0             0.0871             0.1630      0.5032    0.3263   1.0000

# Chl still concurved w Temp0, but all values <0.5

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
#   s(sqrt(CEddyDist0), bs = "cs", k = 4) + s(SSH0, bs = "cs", 
#                                             k = 5) + s(Temp0, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 3.259379   0.002627    1241   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                         edf Ref.df Chi.sq p-value    
#   s(log(abs(FSLE0)))  3.974      4   2412  <2e-16 ***
#   s(log(Chl0))        3.964      4   4217  <2e-16 ***
#   s(Sal700)           3.994      4  13736  <2e-16 ***
#   s(sqrt(CEddyDist0)) 2.997      3   3912  <2e-16 ***
#   s(SSH0)             3.996      4  22785  <2e-16 ***
#   s(Temp0)            3.991      4  15442  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.531   Deviance explained = 58.5%
# -REML =  86053  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: Tweedie(p=1.628) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log_Chl0, bs = "cs", k = 3) + s(Sal400, bs = "cs", k = 5) + 
#   s(sqrt_AEddyDist0, bs = "cs", k = 4) + s(sqrt_CEddyDist0, 
#                                            bs = "cs", k = 3) + s(Temp0, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.90390    0.03534   82.17   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(log_Chl0)        1.951      2 93.34  <2e-16 ***
#   s(Sal400)          3.934      4 31.14  <2e-16 ***
#   s(sqrt_AEddyDist0) 2.913      3 67.29  <2e-16 ***
#   s(sqrt_CEddyDist0) 1.636      2 35.74  <2e-16 ***
#   s(Temp0)           3.654      4 50.55  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.364   Deviance explained =   45%
# -REML = 4280.8  Scale est. = 4.8749    n = 1509

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

