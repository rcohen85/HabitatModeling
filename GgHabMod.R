library(tidyverse)
library(splines2)
library(geepack)
source("getPvalues.R")

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

data = data.frame(read.csv('J:/Chpt_3/ModelData/UD28_masterDF.csv'))
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