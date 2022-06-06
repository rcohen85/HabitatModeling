library(tidyverse)
library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)
library(forecast)
library(nlme)
library(itsadug)

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

## GEEGAM approach ------------------------------------
library(tidyverse)
library(splines2)
library(geepack)
source("getPvalues.R")

data = data.frame(read.csv('J:/Chpt_3/ModelData/UD28_masterDF.csv'))

# # Transform data as necessary
# data$Chl0 = log(data$Chl0)
# data$FSLE0 = log((data$FSLE0^2))
# data$FSLE0[is.infinite(data$FSLE0)] = NA
# 
# # Center & scale covariates
# data[4:dim(data)[2]] = scale(data[4:dim(data)[2]],center=TRUE,scale=TRUE)

# Create grouping variable
lagID = 40
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

# Remove zeros in FSLE data to prepare for later transformation
data$FSLE0[data$FSLE0==0] = NA

# Remove incomplete observations (NAs in FSLE)
badRows = which(is.na(data),arr.ind=TRUE)[,1]
data = data[-badRows,]

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("AEddyDist0",
                  "Chl0",
                  "FSLE0",
                  "Sal0",
                  "Sal200",
                  "EKE0",
                  "SSH0",
                  "Temp0",
                  "Temp200",
                  "VelAsp0",
                  "Slope",
                  "Aspect")

modOpts = c("linMod","threeKnots","fourKnots","fiveKnots")
QIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"VelAsp")){
    periodic = TRUE
  } else { periodic = FALSE}

  modelCall = paste('geeglm(Pres~data$',smoothVarList[i],',data=data,family=poisson,id=GroupID,corstr="ar1")',sep="")
  linMod = eval(parse(text=modelCall))

  modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.5)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')),periodic=FALSE),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  smoothMod1 = eval(parse(text=modelCall))

  modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.333,0.666)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')),periodic=periodic),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  smoothMod2 = eval(parse(text=modelCall))

  modelCall = paste('geeglm(Pres~mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=c(0.275,0.5,0.725)),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')),periodic=periodic),data=data,family=poisson,id=GroupID,corstr=\"ar1\")',sep="")
  smoothMod3 = eval(parse(text=modelCall))

  QIC_votes[i,1:4] = c(QIC(linMod)[[1]],QIC(smoothMod1)[[1]],QIC(smoothMod2)[[1]],QIC(smoothMod3)[[1]])
  QIC_votes[i,5] = modOpts[which.min(QIC_votes[i,1:4])]
}

colnames(QIC_votes) = c(modOpts,"Best")
rownames(QIC_votes) = smoothVarList[]
QIC_votes

#           linMod              threeKnots          fourKnots           Best
# SSH0   "-2517357.26516939" "-2517391.09099838" "-2518046.60409154" "fourKnots"
# Chl0   "-2416276.46810609" "-2460154.55693961" "-2456293.7020221"  "threeKnots"
# Sal0   "-2387459.7791835"  "-2401147.58143417" "-2401909.84094057" "fourKnots"
# Sal200 "-2381759.47681517" "-2398649.87305492" "-2397868.20940183" "threeKnots"
# Temp0  "-2462457.4084357"  "-2503534.24879697" "-2500051.87320411" "threeKnots"
# FSLE0  "-2384991.96877844" "-2385714.26032129" "-2385762.02799539" "fourKnots"
# Slope  "-2409405.0374478"  "-2532732.78179058" "-2533679.36650928" "fourKnots"
# Aspect NA                  NA                  NA                  NA

# Make smooth terms, run full model and check collinearity
smoothVarList = c("SSH0",
                  "Chl0",
                  "Sal0",
                  # "Sal200",
                  "Temp0",
                  # "FSLE0",
                  "Slope")
# "Aspect")
knotList = list(c(0.333,0.666),
                c(0.5),
                c(0.333,0.666),
                # c(0.5),
                c(0.5),
                # c(0.333,0.666),
                c(0.333,0.666))
# c(0.333,0.666))
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
                 corstr="ar1")

VIFvals = vif(fullMod)
VIFvals = cbind(VIFvals,(VIFvals[,3])^2)
colnames(VIFvals)[4] = "LOOK AT ME"
VIFvals

#             GVIF Df GVIF^(1/(2*Df)) LOOK AT ME
# S_SSH0    61.084  5           1.509      2.276
# S_Chl0    11.268  4           1.354      1.832
# S_Sal0    22.547  5           1.366      1.865
# S_Sal200  24.870  4           1.494      2.233
# S_Temp0    8.658  4           1.310      1.715
# S_FSLE0    5.355  5           1.183      1.399
# S_Slope  138.236  5           1.637      2.680
# S_Aspect  13.030  2           1.900      3.610

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

# check term significance
PV = getPvalues(fullMod)

# GETS STUCK THINKING FOREVER (>36HRS)

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

## GAM approach ---------------------
# Regional model
spec = 'SBCD'
outDir = "J:/Chpt_3/GAM_Output"

  # if it doesn't already exist, create directory to save models and figures
  if (!dir.exists(paste(outDir,'/',spec,sep=""))){
    dir.create(paste(outDir,'/',spec,sep=""))
  }

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
  
  summaryData = fullDatesDF %>%
    group_by(WeekID) %>%
    summarize(Pres=sum(Pres,na.rm=TRUE))
  
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

# Remove incomplete observations (NAs in FSLE)
badRows = which(is.na(data),arr.ind=TRUE)[,1]
data = data[-badRows,]
badRows = unique(which(is.na(weeklyDF),arr.ind=TRUE)[,1])
weeklyDF = weeklyDF[-badRows,]


# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("AEddyDist0",
                  "Chl0",
                  "FSLE0",
                  "Sal0",
                  "Sal200",
                  "EKE0",
                  "SSH0",
                  "Temp0",
                  "Temp200",
                  "VelAsp0",
                  "Slope",
                  "Aspect")


# determine whether covars should be included as linear or smooth terms
# seafloor aspect and water velocity direction are included as cyclic smooths

modOpts = c("linMod","threeKnots","fourKnots","fiveKnots")
AIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"Asp")){
    bs = "cc"
  } else { bs = "cs"}
  
  modelCall = paste('gam(Pres~data$',smoothVarList[i],',data=data,family=poisson)',sep="")
  linMod = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=3),data=data,family=poisson)',sep="")
  smoothMod1 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=4),data=data,family=poisson)',sep="")
  smoothMod2 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=5),data=data,family=poisson)',sep="")
  smoothMod3 = eval(parse(text=modelCall))
  
  AIC_votes[i,1:4] = c(AIC(linMod)[[1]],AIC(smoothMod1)[[1]],AIC(smoothMod2)[[1]],AIC(smoothMod3)[[1]])
  AIC_votes[i,5] = modOpts[which.min(AIC_votes[i,1:4])]
}

colnames(AIC_votes) = c(modOpts,"Best")
rownames(AIC_votes) = smoothVarList[]
AIC_votes

#             linMod             threeKnots         fourKnots          fiveKnots          Best       
# AEddyDist0 "476798.02597613"  "473671.83429879"  "473433.202715496" "472644.254851081" "fiveKnots"
# Chl0       "440392.28193109"  "346022.931622249" "306041.418163151" "304705.291256015" "fiveKnots"
# FSLE0      "514061.250033234" "496954.457999153" "496734.39391085"  "494467.442902066" "fiveKnots"
# Sal0       "479826.914056984" "406492.755816828" "382636.328329026" "383238.225164286" "fourKnots"
# Sal200     "469138.179801121" "412433.022344733" "387536.531375634" "387747.705895247" "fourKnots"
# EKE0       "502131.022903349" "488650.358839411" "485664.073084903" "484316.991530492" "fiveKnots"
# SSH0       "311317.800687211" "309939.208716779" "307282.224350419" "305934.601908375" "fiveKnots"
# Temp0      "371767.028384596" "342452.203867807" "341078.002361698" "340585.336411845" "fiveKnots"
# Temp200    "372065.536727268" "350382.520800517" "349370.104565853" "348775.157407256" "fiveKnots"
# VelAsp0    "502304.225654642" "485279.463536039" "485279.463536039" "483942.193107666" "fiveKnots"
# Slope      "465506.280966723" "372437.046930066" "372689.477474422" "370236.922616034" "fiveKnots"
# Aspect     "512688.932044369" "512305.766265523" "512305.766265523" "506574.794585467" "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
              + s(log(Chl0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal200,bs="cs",k=4)
              + s(sqrt(EKE0),bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp200,bs="cs",k=5)
              + s(sqrt(VelAsp0),bs="cc",k=5)
              + s(Slope,bs="cs",k=5)
              + s(sqrt(Aspect),bs="cc",k=5),
              data=data,
              family=poisson,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
              + s(log(Chl0),bs="cs",k=5)
              + s(log(abs(FSLE0)),bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal200,bs="cs",k=4)
              + s(sqrt(EKE0),bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp200,bs="cs",k=5)
              + s(sqrt(VelAsp0),bs="cc",k=5)
              + s(Slope,bs="cs",k=5)
              + s(sqrt(Aspect),bs="cc",k=5),
              data=weeklyDF,
              family=poisson,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

# check convergence
fullMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)
#                       para s(sqrt(AEddyDist0)) s(log(Chl0)) s(log(abs(FSLE0))) s(Sal0) s(Sal200) s(sqrt(EKE0)) s(SSH0) s(Temp0) s(Temp200) s(sqrt(VelAsp0)) s(Slope) s(sqrt(Aspect))
# para                   1              0.0000       0.0000             0.0000  0.0000    0.0000        0.0000  0.0000   0.0000     0.0000           0.0000   0.0000          0.0000
# s(sqrt(AEddyDist0))    0              1.0000       0.0722             0.0575  0.0428    0.0552        0.0402  0.0596   0.0251     0.0636           0.0410   0.0686          0.0928
# s(log(Chl0))           0              0.0728       1.0000             0.0429  0.2874    0.3160        0.0794  0.2121   0.2897     0.3422           0.0787   0.1391          0.0473
# s(log(abs(FSLE0)))     0              0.0508       0.0481             1.0000  0.1525    0.1656        0.0874  0.1212   0.0470     0.0998           0.0855   0.0399          0.0257
# s(Sal0)                0              0.0700       0.2629             0.0938  1.0000    0.9532        0.1915  0.2479   0.2140     0.3581           0.1875   0.1417          0.1014
# s(Sal200)              0              0.0949       0.2619             0.0978  0.9516    1.0000        0.1972  0.2507   0.1990     0.3773           0.1932   0.1741          0.1386
# s(sqrt(EKE0))          0              0.0387       0.0751             0.0724  0.2615    0.2845        1.0000  0.1638   0.0897     0.1637           0.9626   0.0280          0.0467
# s(SSH0)                0              0.1040       0.3703             0.1005  0.4963    0.5295        0.2192  1.0000   0.2985     0.4723           0.2167   0.2114          0.1290
# s(Temp0)               0              0.0462       0.3724             0.0556  0.3451    0.3711        0.1241  0.2669   1.0000     0.7402           0.1218   0.0717          0.0771
# s(Temp200)             0              0.0806       0.3777             0.0932  0.5277    0.5679        0.1823  0.3413   0.7102     1.0000           0.1788   0.0979          0.0798
# s(sqrt(VelAsp0))       0              0.0362       0.0731             0.0713  0.2603    0.2832        0.9222  0.1646   0.0908     0.1637           1.0000   0.0228          0.0471
# s(Slope)               0              0.1038       0.2306             0.0598  0.2204    0.2360        0.0671  0.2838   0.0953     0.1643           0.0647   1.0000          0.2449
# s(sqrt(Aspect))        0              0.1066       0.0552             0.0190  0.1365    0.1360        0.0413  0.1551   0.0379     0.0424           0.0442   0.1831          1.0000

# Chl0 somewhat concurved w SSH, Temp0, Temp200, less w Sal0, Sal200
# Sal0 highly concurved w Sal200, problematic w SSH, Temp0, Temp200, less w Chl0, EKE, VelAsp0
# Sal200 highly concurved w Sal0, prolematic w Chl0, SSH, Temp0, Temp200, EKE, VelAsp0
# EKE problematic w VelAsp0
# Temp200 problematic w Temp0
# Taking out Sal200, VelAsp0, Temp200

conCurv = concurvity(weekMod,full=FALSE)
round(conCurv$estimate,digits=4)
#                       para s(sqrt(AEddyDist0)) s(log(Chl0)) s(log(abs(FSLE0))) s(Sal0) s(Sal200) s(sqrt(EKE0)) s(SSH0) s(Temp0) s(Temp200) s(sqrt(VelAsp0)) s(Slope) s(sqrt(Aspect))
# para                   1              0.0000       0.0000             0.0000  0.0000    0.0000        0.0000  0.0000   0.0000     0.0000           0.0000   0.0000          0.0000
# s(sqrt(AEddyDist0))    0              1.0000       0.0821             0.1184  0.0509    0.0685        0.0435  0.0715   0.0292     0.0751           0.0482   0.0786          0.1047
# s(log(Chl0))           0              0.0811       1.0000             0.0723  0.2852    0.3134        0.0781  0.2120   0.2900     0.3381           0.0857   0.1470          0.0521
# s(log(abs(FSLE0)))     0              0.0849       0.0759             1.0000  0.2117    0.2390        0.1019  0.1837   0.0653     0.1404           0.1183   0.0701          0.0445
# s(Sal0)                0              0.0802       0.2743             0.1478  1.0000    0.9569        0.1622  0.2671   0.2310     0.3613           0.1826   0.1542          0.1120
# s(Sal200)              0              0.1081       0.2720             0.1521  0.9473    1.0000        0.1675  0.2707   0.2097     0.3779           0.1887   0.1944          0.1531
# s(sqrt(EKE0))          0              0.0491       0.1031             0.1427  0.3463    0.3837        1.0000  0.2609   0.1215     0.2120           0.9457   0.0509          0.0976
# s(SSH0)                0              0.1200       0.3717             0.1607  0.5051    0.5454        0.1909  1.0000   0.3026     0.4711           0.2178   0.2217          0.1436
# s(Temp0)               0              0.0493       0.3618             0.0866  0.3654    0.3796        0.1031  0.2631   1.0000     0.7403           0.1199   0.0753          0.0823
# s(Temp200)             0              0.0821       0.3711             0.1485  0.5249    0.5762        0.1527  0.3442   0.7036     1.0000           0.1771   0.1021          0.0861
# s(sqrt(VelAsp0))       0              0.0418       0.1025             0.1341  0.3327    0.3693        0.9263  0.2367   0.1219     0.2117           1.0000   0.0334          0.0677
# s(Slope)               0              0.1239       0.2263             0.1140  0.2098    0.2278        0.0798  0.2892   0.0968     0.1583           0.0789   1.0000          0.2468
# s(sqrt(Aspect))        0              0.1073       0.0547             0.0388  0.1424    0.1432        0.0403  0.1635   0.0402     0.0418           0.0491   0.1846          1.0000
# run reduced model
dayMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
             + s(log(Chl0),bs="cs",k=5)
             + s(log(abs(FSLE0)),bs="cs",k=5)
             + s(Sal0,bs="cs",k=4)
             # + s(Sal200,bs="cs",k=4)
             + s(sqrt(EKE0),bs="cs",k=5)
             + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5),
             # + s(Temp200,bs="cs",k=5)
             # + s(sqrt(VelAsp0),bs="cc",k=5)
             # + s(Slope,bs="cs",k=5)
             # + s(sqrt(Aspect),bs="cc",k=5),
             data=data,
             family=poisson,
             method="REML",
             select=TRUE,
             gamma=1.4,
             na.action="na.fail")

weekMod = gam(Pres ~ s(sqrt(AEddyDist0),bs="cs",k=5)
             + s(log(Chl0),bs="cs",k=5)
             + s(log(abs(FSLE0)),bs="cs",k=5)
             + s(Sal0,bs="cs",k=4)
             # + s(Sal200,bs="cs",k=4)
             + s(sqrt(EKE0),bs="cs",k=5)
             + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5),
             # + s(Temp200,bs="cs",k=5)
             # + s(sqrt(VelAsp0),bs="cc",k=5)
             # + s(Slope,bs="cs",k=5)
             # + s(sqrt(Aspect),bs="cc",k=5),
             data=weeklyDF,
             family=poisson,
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
dayMod$converged
# TRUE
weekMod$converged
# TRUE

# check concurvity of smooth terms
conCurv = concurvity(dayMod,full=FALSE)
round(conCurv$estimate,digits=4)

#                       para s(sqrt(AEddyDist0)) s(log(Chl0)) s(log(abs(FSLE0))) s(Sal0) s(sqrt(EKE0)) s(SSH0) s(Temp0)
# para                   1              0.0000       0.0000             0.0000  0.0000        0.0000  0.0000   0.0000
# s(sqrt(AEddyDist0))    0              1.0000       0.0722             0.0575  0.0428        0.0402  0.0596   0.0251
# s(log(Chl0))           0              0.0728       1.0000             0.0429  0.2874        0.0794  0.2121   0.2897
# s(log(abs(FSLE0)))     0              0.0508       0.0481             1.0000  0.1525        0.0874  0.1212   0.0470
# s(Sal0)                0              0.0700       0.2629             0.0938  1.0000        0.1915  0.2479   0.2140
# s(sqrt(EKE0))          0              0.0387       0.0751             0.0724  0.2615        1.0000  0.1638   0.0897
# s(SSH0)                0              0.1040       0.3703             0.1005  0.4963        0.2192  1.0000   0.2985
# s(Temp0)               0              0.0462       0.3724             0.0556  0.3451        0.1241  0.2669   1.0000

# still some concurvity (esp Sal0 & SSH), but all <0.5

dayModCompTable = dredge(dayMod,
                      beta="none",
                      evaluate=TRUE,
                      trace=TRUE)

weekModCompTable = dredge(weekMod,
                      beta="none",
                      evaluate=TRUE,
                      trace=TRUE)

# run optimal models
optDayMod = get.models(dayModCompTable,subset=1)
optDayMod = optDayMod[[names(optDayMod)]]
save(optDayMod,dayModCompTable,file=paste(outDir,'/',spec,'/','DailyRegionalModel.Rdata',sep=""))
optWeekMod = get.models(weekModCompTable,subset=1)
optWeekMod = optWeekMod[[names(optWeekMod)]]
save(optWeekMod,weekModCompTable,file=paste(outDir,'/',spec,'/','WeeklyRegionalModel.Rdata',sep=""))

# check p-values
summary(optDayMod)

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
# (Intercept) 5.194514   0.002625    1979   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(log(abs(FSLE0)))  3.983      4   2164  <2e-16 ***
#   s(log(Chl0))        3.975      4   3979  <2e-16 ***
#   s(Sal0)             2.982      3  10316  <2e-16 ***
#   s(sqrt(AEddyDist0)) 3.956      4   3272  <2e-16 ***
#   s(sqrt(EKE0))       3.925      4   1053  <2e-16 ***
#   s(SSH0)             3.994      4  15369  <2e-16 ***
#   s(Temp0)            3.989      4   8314  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.709   Deviance explained = 75.6%
# -REML =  36447  Scale est. = 1         n = 1509

# plot
png(filename=paste(outDir,'/',spec,'/',spec,'_allSitesDaily.png',sep=""),width=600,height=600)
plot.gam(optDayMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
while (dev.cur()>1) {dev.off()}
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
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'.png',sep=""),width=600,height=600)
    plot.gam(siteDayModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
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
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'.png',sep=""),width=600,height=600)
    plot.gam(siteWeekModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
    while (dev.cur()>1) {dev.off()}
  }
}

save(siteDayModList,pValDayList,siteDayModCompList,file=paste(outDir,'/',spec,'/','DailySiteSpecificModels.Rdata',sep=""))
save(siteWeekModList,pValWeekList,siteWeekModCompList,file=paste(outDir,'/',spec,'/','WeeklySiteSpecificModels.Rdata',sep=""))

