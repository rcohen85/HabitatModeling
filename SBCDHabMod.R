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

# library(lubridate)
# library(car)
# library(splines2)
# library(ggplot2)
# library(gridExtra)
# library(zoo)
# library(pracma)
# library(mgcv)
# library(SimDesign)
# library(multitaper)

spec = 'SBCD'
data = data.frame(read.csv('J:/Chpt_3/ModelData/UD28_masterDF.csv'))
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
smoothVarList = c("sqrt_AEddyDist0",
                  "sqrt_CEddyDist0",
                  "log_Chl0",
                  "Sal0",
                  "Sal200",
                  "SSH0",
                  "Temp0",
                  "Temp200")

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
                   + bs(log_Chl0)
                   + bs(Sal0)
                   + bs(Sal200)
                   + bs(SSH0)
                   + bs(Temp0)
                   + bs(Temp200),
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
# WC     3
# NFC    3
# HAT    3
# GS     2
# BP     2
# BS     2

# test for overdispersion
dispMod = glm(Pres~bs(sqrt_AEddyDist0) # include all terms in smoothVarList above!!
              + bs(sqrt_CEddyDist0)
              + bs(log_Chl0)
              + bs(Sal0)
              + bs(Sal200)
              + bs(SSH0)
              + bs(Temp0)
              + bs(Temp200),
              data=weeklyDF,family=poisson)

dispersiontest(dispMod,alternative='two.sided')

# Dispersion test
# 
# data:  dispMod
# z = 17.724, p-value < 2.2e-16
# alternative hypothesis: true dispersion is not equal to 1
# sample estimates:
#   dispersion 
# 60.32127

# data are very overdispersed, will use Tweedie family in models
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

#                 linMod             threeKnots         fourKnots          fiveKnots          Best       
# sqrt_AEddyDist0 "19976.5736094526" "19976.7023662434" "19969.6862227397" "19971.3922233781" "fourKnots"
# sqrt_CEddyDist0 "20084.188031296"  "20084.3208779169" "20084.4735173124" "20084.5928875426" "linMod"   
# log_Chl0        "18956.257003049"  "18777.5374430764" "18668.5029501613" "18645.7403004278" "fiveKnots"
# Sal0            "19959.4187959122" "19564.4908709053" "19420.0636529896" "19425.4680005534" "fourKnots"
# Sal200          "19905.4621213139" "19605.7179791146" "19439.2560028503" "19440.7616918607" "fourKnots"
# SSH0            "18698.3346469256" "18684.5517893375" "18684.0440361589" "18616.3983001702" "fiveKnots"
# Temp0           "19421.5824478787" "19267.7805932054" "19256.1727988362" "19254.6535885851" "fiveKnots"
# Temp200         "19393.8130514412" "19281.4796389458" "19277.7841525424" "19276.9231939763" "fiveKnots"

# run full model
# fullMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=4)
#               + s(log(Chl0),bs="cs",k=5)
#               + s(log(abs(FSLE0)),bs="cs",k=5)
#               + s(Sal0,bs="cs",k=4)
#               + s(Sal200,bs="cs",k=4)
#               + s(sqrt(EKE0),bs="cs",k=5)
#               + s(SSH0,bs="cs",k=5)
#               + s(Temp0,bs="cs",k=5)
#               + s(Temp200,bs="cs",k=5)
#               + s(sqrt(VelAsp0),bs="cc",k=5),
#               data=data,
#               family=poisson,
#               method="REML",
#               select=TRUE,
#               gamma=1.4,
#               na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt_AEddyDist0,bs="cs",k=4)
              + sqrt_CEddyDist0
              + s(log_Chl0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal200,bs="cs",k=4)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp200,bs="cs",k=5),
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
#                     para s(sqrt_AEddyDist0) s(log_Chl0) s(Sal0) s(Sal200) s(SSH0) s(Temp0) s(Temp200)
# para                  1             0.0000      0.0000  0.0000    0.0000  0.0000   0.0000     0.0000
# s(sqrt_AEddyDist0)    0             1.0000      0.0804  0.0490    0.0664  0.0626   0.0284     0.0733
# s(log_Chl0)           0             0.1502      1.0000  0.2852    0.3134  0.2120   0.2900     0.3381
# s(Sal0)               0             0.1486      0.2743  1.0000    0.9569  0.2671   0.2310     0.3613
# s(Sal200)             0             0.2032      0.2720  0.9473    1.0000  0.2707   0.2097     0.3779
# s(SSH0)               0             0.1850      0.3717  0.5051    0.5454  1.0000   0.3026     0.4711
# s(Temp0)              0             0.0871      0.3618  0.3654    0.3796  0.2631   1.0000     0.7403
# s(Temp200)            0             0.1506      0.3711  0.5249    0.5762  0.3442   0.7036     1.0000

# Sal0 highly concurved w Sal200, problematic w SSH, Temp200
# Sal200 highly concurved w Sal0, problematic w SSH, Temp200
# Temp200 problematic w Temp0

# Taking out Sal0,  Temp200

# dayMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
#              + s(log(Chl0),bs="cs",k=5)
#              + s(log(abs(FSLE0)),bs="cs",k=5)
#              + s(Sal0,bs="cs",k=4)
#              # + s(Sal200,bs="cs",k=4)
#              + s(sqrt(EKE0),bs="cs",k=5)
#              + s(SSH0,bs="cs",k=5)
#              + s(Temp0,bs="cs",k=5),
#              # + s(Temp200,bs="cs",k=5)
#              # + s(sqrt(VelAsp0),bs="cc",k=5),
#              data=data,
#              family=poisson,
#              method="REML",
#              select=TRUE,
#              gamma=1.4,
#              na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt_AEddyDist0,bs="cs",k=4)
              + sqrt_CEddyDist0
              + s(log_Chl0,bs="cs",k=5)
              # + s(Sal0,bs="cs",k=4)
              + s(Sal200,bs="cs",k=4)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5),
              # + s(Temp200,bs="cs",k=5),
             data=weeklyDF,
             family=modFam,
             method="REML",
             select=TRUE,
             gamma=1.4,
             na.action="na.fail")


# # Model autcorrelation of presence time series for inclusion in models
# fit= auto.arima(data$Pres,seasonal=TRUE)
# p <- length(fit$model$phi)
# q <- length(fit$model$theta)
# 
# # specify autocorrelated chunks (sites)
# data = start_event(data,column="Date",event="Site")
# 
# # determine autcorrelation value of 1st lag
# r1 <- start_value_rho(fullMod, plot=TRUE)
# 
# corMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
#              + s(log(Chl0),bs="cs",k=5)
#              + s(log(abs(FSLE0)),bs="cs",k=5)
#              + s(Sal0,bs="cs",k=4)
#              # + s(Sal200,bs="cs",k=4)
#              + s(sqrt(EKE0),bs="cs",k=5)
#              + s(SSH0,bs="cs",k=5)
#              + s(Temp0,bs="cs",k=5)
#              # + s(Temp200,bs="cs",k=5)
#              # + s(sqrt(VelAsp0),bs="cc",k=5)
#              + s(Slope,bs="cs",k=5)
#              + s(sqrt(Aspect),bs="cc",k=5),
#              data=data,
#              family=poisson,
#              corstr=corARMA(p=p,q=q),
#              method="REML",
#              select=TRUE,
#              gamma=1.4,
#              na.action="na.fail")
# 
# corMod2 = bam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
#              + s(log(Chl0),bs="cs",k=5)
#              + s(log(abs(FSLE0)),bs="cs",k=5)
#              + s(Sal0,bs="cs",k=4)
#              # + s(Sal200,bs="cs",k=4)
#              + s(sqrt(EKE0),bs="cs",k=5)
#              + s(SSH0,bs="cs",k=5)
#              + s(Temp0,bs="cs",k=5)
#              # + s(Temp200,bs="cs",k=5)
#              # + s(sqrt(VelAsp0),bs="cc",k=5)
#              + s(Slope,bs="cs",k=5)
#              + s(sqrt(Aspect),bs="cc",k=5),
#              data=data,
#              family=poisson,
#              rho=r1,
#              AR.start=data$start.event,
#              method="REML",
#              select=TRUE,
#              gamma=1.4,
#              na.action="na.fail")

# check convergence
# dayMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)
#                     para s(sqrt_AEddyDist0) s(log_Chl0) s(Sal200) s(SSH0) s(Temp0)
# para                  1             0.0000      0.0000    0.0000  0.0000   0.0000
# s(sqrt_AEddyDist0)    0             1.0000      0.0804    0.0664  0.0626   0.0284
# s(log_Chl0)           0             0.1502      1.0000    0.3134  0.2120   0.2900
# s(Sal200)             0             0.2032      0.2720    1.0000  0.2707   0.2097
# s(SSH0)               0             0.1850      0.3717    0.5454  1.0000   0.3026
# s(Temp0)              0             0.0871      0.3618    0.3796  0.2631   1.0000

# still some concurvity btwn Sal200 and SSH, but proceeding

# dayModCompTable = dredge(dayMod,
#                       beta="none",
#                       evaluate=TRUE,
#                       trace=TRUE)

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
#   Pres ~ s(log(abs(FSLE0)), bs = "cs", k = 5) + s(log(Chl0), bs = "cs", 
#                                                   k = 5) + s(Sal0, bs = "cs", k = 4) + s(sqrt(AEddyDist0), 
#                                                                                          bs = "cs", k = 5) + s(sqrt(EKE0), bs = "cs", k = 5) + s(SSH0, 
#                                                                                                                                                  bs = "cs", k = 5) + s(Temp0, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 3.264665   0.002614    1249   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq p-value    
# s(log(abs(FSLE0)))  3.964      4  2356.8  <2e-16 ***
#   s(log(Chl0))        3.969      4  4969.2  <2e-16 ***
#   s(Sal0)             2.988      3 13913.3  <2e-16 ***
#   s(sqrt(AEddyDist0)) 3.973      4  3436.5  <2e-16 ***
#   s(sqrt(EKE0))       3.807      4   248.4  <2e-16 ***
#   s(SSH0)             3.994      4 18627.7  <2e-16 ***
#   s(Temp0)            3.991      4 10553.5  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.536   Deviance explained = 58.7%
# -REML =  85804  Scale est. = 1         n = 10375

summary(optWeekMod)

# Family: Tweedie(p=1.53) 
# Link function: log 
# 
# Formula:
#   Pres ~ s(log_Chl0, bs = "cs", k = 5) + s(Sal200, bs = "cs", k = 4) + 
#   s(sqrt_AEddyDist0, bs = "cs", k = 4) + s(SSH0, bs = "cs", 
#                                            k = 5) + s(Temp0, bs = "cs", k = 5) + sqrt_CEddyDist0 + 1
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      5.383357   0.073624  73.120   <2e-16 ***
#   sqrt_CEddyDist0 -0.014956   0.005967  -2.506   0.0123 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df       F p-value    
# s(log_Chl0)        3.248      4  18.663  <2e-16 ***
#   s(Sal200)          2.689      3  96.767  <2e-16 ***
#   s(sqrt_AEddyDist0) 1.289      3   9.265  <2e-16 ***
#   s(SSH0)            3.776      4 109.331  <2e-16 ***
#   s(Temp0)           2.856      4  34.699  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.697   Deviance explained = 73.5%
# -REML = 6438.7  Scale est. = 3.5865    n = 1509

# plot
# png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesDaily.png',sep=""),width=600,height=600)
# plot.gam(optDayMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0,residuals=TRUE)
# while (dev.cur()>1) {dev.off()}
png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesWeekly.png',sep=""),width=600,height=600)
plot.gam(optWeekMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0,residuals=FALSE)
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
    fullSiteDayMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
                         + s(log(Chl0),bs="cs",k=5)
                         + s(log(abs(FSLE0)),bs="cs",k=5)
                         + s(Sal0,bs="cs",k=4)
                         + s(sqrt(EKE0),bs="cs",k=5)
                         + s(SSH0,bs="cs",k=5)
                         + s(Temp0,bs="cs",k=5),
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
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'_DailyR.png',sep=""),width=600,height=600)
    plot.gam(siteDayModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0,residuals=TRUE)
    while (dev.cur()>1) {dev.off()}
    
  }
  
  if (sum(weekData$Pres>0)>25){
    #run full weekly model for this site
    fullSiteWeekMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
                          + s(log(Chl0),bs="cs",k=5)
                          + s(log(abs(FSLE0)),bs="cs",k=5)
                          + s(Sal0,bs="cs",k=4)
                          + s(sqrt(EKE0),bs="cs",k=5)
                          + s(SSH0,bs="cs",k=5)
                          + s(Temp0,bs="cs",k=5),
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
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'_WeeklyR.png',sep=""),width=600,height=600)
    plot.gam(siteWeekModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0,residuals=TRUE)
    while (dev.cur()>1) {dev.off()}
  }
}

save(siteDayModList,pValDayList,siteDayModCompList,file=paste(outDir,'/',spec,'/','DailySiteSpecificModels.Rdata',sep=""))
save(siteWeekModList,pValWeekList,siteWeekModCompList,file=paste(outDir,'/',spec,'/','WeeklySiteSpecificModels.Rdata',sep=""))



## GEEGAM approach ------------------------------------
library(tidyverse)
library(splines2)
library(splines)
library(geepack)
library(ggfortify)
library(car)
source("getPvalues.R")

outDir = "J:/Chpt_3/GEEGLM_Output"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

# determine block size based on autocorrelation
sites = unique(weeklyDF$Site)
siteCorr = data.frame(Site=character(length=length(sites)),Corr=numeric(length=length(sites)))
for (i in 1:length(sites)){
  
  siteInd = which(weeklyDF$Site==sites[i])
  BlockMod = glm(Pres~bs(Temp0)
                 + bs(Sal0)
                 + bs(log_Chl0)
                 + bs(log_abs_FSLE0)
                 + bs(VelMag0)
                 + bs(sqrt_VelAsp0)
                 + bs(sqrt_EKE0)
                 + bs(SSH0)
                 + bs(sqrt_AEddyDist0)
                 + bs(sqrt_CEddyDist0)
                 + bs(GSDist_div100),
                 data=weeklyDF[siteInd,],family=poisson)
  
  acorr = acf(residuals(BlockMod), lag.max = 100, main=paste(spec,"at",sites[i]))
  CI = ggfortify:::confint.acf(acorr)
  ACFidx = which(acorr[["acf"]] < CI, arr.ind=TRUE)
  siteCorr$Site[i] = sites[i]
  siteCorr$Corr[i] = ACFidx[1]
}

# Create grouping variable
lagID = max(siteCorr$Corr)
numClust = length(weeklyDF$Pres)/(lagID-1)
if (numClust<length(weeklyDF$Pres)){
  clustID = rep(1:ceiling(numClust),each=lagID)
  clustID = clustID[1:length(weeklyDF$Pres)]
} else {
  clustID = 1:length(weeklyDF$Pres)
}
weeklyDF$GroupID = clustID

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("sqrt_AEddyDist0",
                  "log_Chl0",
                  "log_abs_FSLE0",
                  "Sal0",
                  "Sal200",
                  "sqrt_EKE0",
                  "SSH0",
                  "Temp0",
                  "Temp200",
                  "sqrt_VelAsp0")

modOpts = c("linMod","threeKnots","fourKnots","fiveKnots")
QIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"VelAsp")){
    periodic = TRUE
  } else { periodic = FALSE}
  
  modelCall = paste('geeglm(Pres~weeklyDF$',smoothVarList[i],',data=weeklyDF,family=poisson,id=GroupID,corstr="ar1")',sep="")
  linMod = eval(parse(text=modelCall))
  
  modelCall = paste('geeglm(Pres~mSpline(weeklyDF$',smoothVarList[i],',knots=quantile(weeklyDF$',smoothVarList[i],',probs=c(0.5)),Boundary.knots=c(min(weeklyDF$',smoothVarList[i],'),max(weeklyDF$',smoothVarList[i],')),periodic=FALSE),data=weeklyDF,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  smoothMod1 = eval(parse(text=modelCall))
  
  modelCall = paste('geeglm(Pres~mSpline(weeklyDF$',smoothVarList[i],',knots=quantile(weeklyDF$',smoothVarList[i],',probs=c(0.333,0.666)),Boundary.knots=c(min(weeklyDF$',smoothVarList[i],'),max(weeklyDF$',smoothVarList[i],')),periodic=periodic),data=weeklyDF,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  smoothMod2 = eval(parse(text=modelCall))
  
  modelCall = paste('geeglm(Pres~mSpline(weeklyDF$',smoothVarList[i],',knots=quantile(weeklyDF$',smoothVarList[i],',probs=c(0.275,0.5,0.725)),Boundary.knots=c(min(weeklyDF$',smoothVarList[i],'),max(weeklyDF$',smoothVarList[i],')),periodic=periodic),data=weeklyDF,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  smoothMod3 = eval(parse(text=modelCall))
  
  QIC_votes[i,1:4] = c(QIC(linMod)[[1]],QIC(smoothMod1)[[1]],QIC(smoothMod2)[[1]],QIC(smoothMod3)[[1]])
  QIC_votes[i,5] = modOpts[which.min(QIC_votes[i,1:4])]
}

colnames(QIC_votes) = c(modOpts,"Best")
rownames(QIC_votes) = smoothVarList[]
QIC_votes

#                   linMod              threeKnots          fourKnots           fiveKnots           Best        
# sqrt_AEddyDist0 "-4077848.80747756" "-4069019.4742813"  "-4071331.43807784" "-4070865.15796378" "linMod"    
# log_Chl0        "-4197691.23175159" "-4242978.40861398" "-4245436.890391"   "-4244611.93511617" "fourKnots" 
# log_abs_FSLE0   "-4071612.81747493" "-4073962.42515323" "-4073898.53231251" "-4073847.43977969" "threeKnots"
# Sal0            "-4099184.76027428" "-4160037.25766481" "-4163786.89124787" "-4162471.88363391" "fourKnots" 
# Sal200          "-4100360.37337247" "-4163212.29616996" "-4167107.78406302" "-4165991.1842822"  "fourKnots" 
# sqrt_EKE0       "-4072061.30467109" "-4076196.37131878" "-4076292.80972148" "-4076519.89228992" "fiveKnots" 
# SSH0            "-4269029.57039468" "-4273980.18594461" "-4280732.43017986" "-4277099.39312909" "fourKnots" 
# Temp0           "-4202891.2950828"  "-4237853.79771032" "-4235690.00177933" "-4235527.42574321" "threeKnots"
# Temp200         "-4202466.8722529"  "-4235564.42989713" "-4235931.78258385" "-4236115.14178248" "fiveKnots" 
# sqrt_VelAsp0    "-4072515.86589221" "-4078138.1611877"  "-4078023.43585877" "-4078133.24259284" "threeKnots"

# Make smooth terms, run full model and check collinearity
smoothVarList = c("log_Chl0",
                  "log_abs_FSLE0",
                  "Sal0",
                  "Sal200",
                  "sqrt_EKE0",
                  "SSH0",
                  "Temp0",
                  "Temp200",
                  "sqrt_VelAsp0")

knotList = list(c(0.333,0.666),
                c(0.5),
                c(0.333,0.666),
                c(0.333,0.666),
                c(0.275,0.5,0.725),
                c(0.333,0.666),
                c(0.5),
                c(0.275,0.5,0.725),
                c(0.333,0.666)) # velocity aspect has to have at least 4 knots to be circular

linVarList = list("sqrt_AEddyDist0")
smoothNameList = character()


for (i in 1:length(smoothVarList)){
  
  if (str_detect(smoothVarList[i],"Asp")){ periodic=TRUE} else {periodic=FALSE}
  eval(parse(text=paste('S_',smoothVarList[i],'= mSpline(weeklyDF$',smoothVarList[i],',knots=quantile(weeklyDF$',smoothVarList[i],',probs=unlist(knotList[i])),Boundary.knots=c(min(weeklyDF$',smoothVarList[i],'),max(weeklyDF$',smoothVarList[i],')),periodic=periodic)',sep="")))
  
  smoothNameList = c(smoothNameList,paste('S_',smoothVarList[i],sep=""))
}
thisForm = formula(paste('Pres~',paste(c(smoothNameList,linVarList),collapse="+"),sep=""))

fullMod = geeglm(thisForm,
                 family=poisson,
                 data=weeklyDF,
                 id=GroupID,
                 corstr="ar1",
                 na.action="na.fail")

# check convergence
fullMod$geese$error
# 0 model converged

# check collinearity
VIFvals = vif(fullMod)
VIFvals = cbind(VIFvals,(VIFvals[,3])^2)
colnames(VIFvals)[4] = "LOOK AT ME"
VIFvals

#                         GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# S_log_Chl0      1.565263e+01  5        1.316615   1.733474
# S_log_abs_FSLE0 6.781570e+00  4        1.270329   1.613736
# S_Sal0          2.431228e+07  5        5.477502  30.003033
# S_Sal200        1.140465e+07  5        5.078181  25.787927
# S_sqrt_EKE0     6.541867e+02  6        1.716492   2.946345
# S_SSH0          9.929238e+01  5        1.583768   2.508321
# S_Temp0         5.118689e+03  4        2.908338   8.458429
# S_Temp200       3.250016e+04  6        2.376788   5.649121
# S_sqrt_VelAsp0  1.321150e+02  2        3.390299  11.494130
# sqrt_AEddyDist0 1.651042e+00  1        1.284929   1.651042

# removing Sal0
fullMod = update(fullMod,. ~ . - S_Sal0)

#                         GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# S_log_Chl0        14.471109  5        1.306322   1.706477
# S_log_abs_FSLE0    4.565246  4        1.209019   1.461726
# S_Sal200          97.572012  5        1.581002   2.499569
# S_sqrt_EKE0      580.814858  6        1.699560   2.888504
# S_SSH0            88.904848  5        1.566363   2.453495
# S_Temp0         1738.467508  4        2.541095   6.457161
# S_Temp200       9671.593774  6        2.148448   4.615829
# S_sqrt_VelAsp0   123.121508  2        3.331068  11.096013
# sqrt_AEddyDist0    1.392757  1        1.180151   1.392757

# removing VelAsp0
fullMod = update(fullMod,. ~ . - S_sqrt_VelAsp0)

#                         GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# S_log_Chl0        13.899573  5        1.301068   1.692779
# S_log_abs_FSLE0    4.057004  4        1.191312   1.419225
# S_Sal200          91.486601  5        1.570854   2.467581
# S_sqrt_EKE0        7.461110  6        1.182316   1.397871
# S_SSH0            84.723268  5        1.558835   2.429968
# S_Temp0         1372.938782  4        2.467212   6.087133
# S_Temp200       7648.395135  6        2.106837   4.438761
# sqrt_AEddyDist0    1.385384  1        1.177023   1.385384

# removing Temp200
fullMod = update(fullMod,. ~ . - S_Temp200)

#                       GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# S_log_Chl0      11.680513  5        1.278634   1.634904
# S_log_abs_FSLE0  3.103345  4        1.152070   1.327265
# S_Sal200        32.550266  5        1.416627   2.006831
# S_sqrt_EKE0      4.937486  6        1.142332   1.304921
# S_SSH0          67.556675  5        1.523936   2.322382
# S_Temp0          9.379330  4        1.322883   1.750020
# sqrt_AEddyDist0  1.330831  1        1.153616   1.330831

# check convergence
fullMod$geese$error
# 0 model converged

# dredge model
GEECompTable = dredge(fullMod,
                      beta="none",
                      evaluate=TRUE,
                      trace=TRUE)

# get optimal model
optMod = get.models(GEECompTable,subset=1)
optMod = optMod[[names(optMod)]]
save(optMod,GEECompTable,file=paste(outDir,'/',spec,'/','WeeklyRegionalModel.Rdata',sep=""))

# # check term significance
# PV = getPvalues(optMod)
# #           Variable  p-value
# # 1      S_log_Chl0  <0.0001
# # 2 S_log_abs_FSLE0  0.02229
# # 3        S_Sal200  <0.0001
# # 4     S_sqrt_EKE0 0.447738
# # 5          S_SSH0  <0.0001
# # 6         S_Temp0  <0.0001
# # 7 sqrt_AEddyDist0 0.129326
# 
# #Remove non-significant terms, re-run model, check p-values
# PV$'p-value'[PV$'p-value'=="<0.0001"] = 0.0001
# badVars = PV$Variable[as.numeric(PV$'p-value')>=0.05]
# redMod<-eval(parse(text=paste("update(fullMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
# if (redMod$geese$error==1){
#   print("Model did not converge")
# } else {PVred = getPvalues(redMod)
# PVred$'p-value'[PVred$'p-value'=="<0.0001"] = 0.0001}
# 
# #           Variable  p-value
# # 1      S_log_Chl0  <0.0001
# # 2 S_log_abs_FSLE0 0.037413
# # 3        S_Sal200 0.000203
# # 4          S_SSH0  <0.0001
# # 5         S_Temp0  <0.0001

# Plot terms from regional model
source("plotSmooths.R")
source("plotLinears.R")
terms = names(optMod$model)[2:length(names(optMod$model))]
plotList = list()
resids=TRUE

for (i in 1:length(terms)){
  if (str_detect(terms[i],"S_")){ # plot smooth terms
    term = str_remove(terms[i],"S_")
    knotInd = which(!is.na(str_match(smoothVarList,term)))
    k=length(unlist(knotList[knotInd]))
    coefInd = which(str_detect(names(optMod$coefficients),term))
    if (str_detect(term,"Asp")){periodic=TRUE} else {periodic=FALSE}
    plotList[[term]] = print(plotSmooths(optMod,weeklyDF,term,coefInd,k,periodic,resids,site=NA,title=NULL))
  } else { # plot linear terms
    term=terms[i]
    coefInd = which(str_detect(names(optMod$coefficients),term))
    plotList[[term]] = print(plotLinears(optMod,weeklyDF,term,coefInd,resids,site=NA,title=NULL))
  }
}
png(file=paste(outDir,'/',spec,'/',spec,'_allSitesWeeklyRlog10.png',sep=""),width=800,height=700,)
grid.arrange(grobs=plotList,nrow=3)
while (dev.cur()>1) {dev.off()}

# site-specific models ---------------------------------------------------
siteWeekModList = list()
pValWeekList = list()
siteWeekModCompList = list()
thisForm = formula(fullMod)
resids = TRUE

for (i in 1:length(sites)){
  
  weekInd = which(!is.na(str_match(weeklyDF$Site,sites[i])))
  weekData = weeklyDF[weekInd,]
  
  if (sum(weekData$Pres>0)>25){
    
    # remake smooth terms with data just from this site
    for (k in 1:length(smoothVarList)){
      if (str_detect(smoothVarList[k],"Asp")){ periodic=TRUE} else {periodic=FALSE}
      eval(parse(text=paste('S_',smoothVarList[k],'= mSpline(weeklyDF$',smoothVarList[k],'[weekInd],knots=quantile(weeklyDF$',smoothVarList[k],'[weekInd],probs=unlist(knotList[k])),Boundary.knots=c(min(weeklyDF$',smoothVarList[k],'[weekInd]),max(weeklyDF$',smoothVarList[k],'[weekInd])),periodic=periodic)',sep="")))
    }
    
    #run full weekly model for this site
    fullSiteWeekMod = geeglm(thisForm,
                             family=poisson,
                             data=weekData,
                             id=GroupID,
                             corstr="ar1",
                             na.action="na.fail")
    
    siteWeekModCompTable = dredge(fullSiteWeekMod,
                                  beta="none",
                                  evaluate=TRUE,
                                  trace=TRUE)
    siteWeekModCompList[[sites[i]]] = siteWeekModCompTable
    
    optSiteWeekMod = get.models(siteWeekModCompTable,subset=1)
    optSiteWeekMod = optSiteWeekMod[[names(optSiteWeekMod)]]
    siteWeekPV = summary(optSiteWeekMod)$s.pv
    
    siteWeekModList[[sites[i]]] = optSiteWeekMod
    pValWeekList[[sites[i]]] = siteWeekPV
    
    sink(paste(outDir,'/',spec,'/',spec,'_',sites[i],'_WeeklySummary.txt',sep=""))
    print(summary(siteWeekModList[[sites[i]]]))
    sink()
    
    terms = names(optSiteWeekMod$model)[2:length(names(optSiteWeekMod$model))]
    plotList = list()
    
    for (j in 1:length(terms)){
      if (str_detect(terms[j],"S_")){ # plot smooth terms
        term = str_remove(terms[j],"S_")
        knotInd = which(!is.na(str_match(smoothVarList,term)))
        k=length(unlist(knotList[knotInd]))
        coefInd = which(str_detect(names(optSiteWeekMod$coefficients),term))
        if (str_detect(term,"Asp")){periodic=TRUE} else {periodic=FALSE}
        plotList[[term]] = print(plotSmooths(optSiteWeekMod,weekData,term,coefInd,k,periodic,resids,site=NA,title=NULL))
      } else { # plot linear terms
        term=terms[j]
        coefInd = which(str_detect(names(optSiteWeekMod$coefficients),term))
        plotList[[term]] = print(plotLinears(optSiteWeekMod,weekData,term,coefInd,resids,site=NA,title=NULL))
      }
    }
    png(file=paste(outDir,'/',spec,'/',spec,'_',sites[i],'_Weekly.png',sep=""),width=800,height=700,)
    grid.arrange(grobs=plotList,nrow=3)
    while (dev.cur()>1) {dev.off()}
  }
  
}

save(siteWeekModList,pValWeekList,siteWeekModCompList,file=paste(outDir,'/',spec,'/','WeeklySiteSpecificGEEs.Rdata',sep=""))

