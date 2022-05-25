library(tidyverse)
# library(car)
library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)


# library(tidyverse)
# library(splines2)
# library(geepack)
# source("getPvalues.R")

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
data = data.frame(read.csv('J:/Chpt_3/ModelData/UD28_masterDF.csv'))
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

# Test for how a term should be included in the model
startTime = Sys.time()
smoothVarList = c("SSH0",
                  "Chl0",
                  "Sal0",
                  "Sal200",
                  "Temp0",
                  "FSLE0",
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

data = data.frame(read.csv('J:/Chpt_3/ModelData/UD28_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

# Check which covars are correlated w presence to determine starting covar list

# Test for how a term should be included in the model
smoothVarList = c("Chl0",
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
AIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"Slope")|str_detect(smoothVarList[i],"Asp")){
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
# Chl0    "446062.513080801" "351317.638746874" "310334.178303824" "309050.026050442" "fiveKnots"
# FSLE0   "520088.534596979" "502724.700277915" "502631.878618275" "500301.060225707" "fiveKnots"
# Sal0    "486235.085207736" "411360.941415359" "387124.722844761" "387743.671350695" "fourKnots"
# Sal200  "475574.470490484" "416908.505680563" "392087.760434819" "392233.025805463" "fourKnots"
# EKE0    "508364.452598272" "494830.704327545" "491753.828970785" "490354.459864932" "fiveKnots"
# SSH0    "315223.54295346"  "313716.500031201" "310937.15943765"  "309552.283203949" "fiveKnots"
# Temp0   "376639.997067505" "346409.991686607" "345047.680730336" "344559.370920586" "fiveKnots"
# Temp200 "377232.102712691" "354782.226873767" "353815.24837876"  "353224.407000695" "fiveKnots"
# VelAsp0 "508540.733795566" "491357.274370519" "491357.274370519" "490008.048378874" "fiveKnots"
# Slope   "471243.092863463" "378355.828929701" "378355.828929701" "378321.313830591" "fiveKnots"
# Aspect  "519177.311573147" "518337.881472128" "518337.881472128" "512187.826366409" "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
              + s(FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal200,bs="cs",k=4)
              + s(EKE0,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp200,bs="cs",k=5)
              + s(VelAsp0,bs="cc",k=5)
              + s(Slope,bs="cc",k=5)
              + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail")

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)

#             para s(Chl0) s(FSLE0) s(Sal0) s(Sal200) s(EKE0) s(SSH0) s(Temp0) s(Temp200) s(VelAsp0) s(Slope) s(Aspect)
# para          1  0.0000   0.0000  0.0000    0.0000  0.0000  0.0000   0.0000     0.0000     0.0000   0.0000    0.0000
# s(Chl0)       0  1.0000   0.0322  0.3345    0.3936  0.0776  0.2111   0.2806     0.3368     0.0699   0.1455    0.0477
# s(FSLE0)      0  0.0206   1.0000  0.1710    0.1938  0.0812  0.1191   0.0448     0.0974     0.0712   0.0450    0.0267
# s(Sal0)       0  0.1286   0.1268  1.0000    0.9617  0.1797  0.2476   0.2139     0.3567     0.1608   0.1521    0.1053
# s(Sal200)     0  0.1190   0.1324  0.9562    1.0000  0.1852  0.2503   0.1988     0.3761     0.1657   0.1784    0.1463
# s(EKE0)       0  0.0315   0.1023  0.2900    0.3294  1.0000  0.1574   0.0864     0.1597     0.9430   0.0335    0.0395
# s(SSH0)       0  0.1551   0.1434  0.5754    0.6569  0.2027  1.0000   0.2993     0.4723     0.1849   0.2145    0.1298
# s(Temp0)      0  0.2332   0.0840  0.3741    0.4189  0.1160  0.2677   1.0000     0.7417     0.1063   0.0784    0.0809
# s(Temp200)    0  0.2029   0.1263  0.5911    0.6660  0.1699  0.3421   0.7125     1.0000     0.1537   0.1090    0.0840
# s(VelAsp0)    0  0.0318   0.0981  0.2848    0.3241  0.8862  0.1511   0.0873     0.1588     1.0000   0.0298    0.0366
# s(Slope)      0  0.0856   0.0484  0.2314    0.2731  0.0551  0.1906   0.0881     0.1557     0.0508   1.0000    0.1961
# s(Aspect)     0  0.0361   0.0194  0.1489    0.1454  0.0348  0.1523   0.0382     0.0431     0.0419   0.1538    1.0000

# Sal0 problematic w SSH, Sal200, Temp0, less w Chl0
# Sal200 prolematic w Chl0, Sal0, SSH, Temp0, Temp200
# EKE problematic w VelAsp0
# Temp200 problematic w Temp0
# Taking out Sal200, VelAsp0, Temp200

# run reduced model
redMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
             + s(FSLE0,bs="cs",k=5)
             + s(Sal0,bs="cs",k=4)
             # + s(Sal200,bs="cs",k=4)
             + s(EKE0,bs="cs",k=5)
             + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5)
             # + s(Temp200,bs="cs",k=5)
             # + s(VelAsp0,bs="cc",k=5)
             + s(Slope,bs="cc",k=5)
             + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail")

# check concurvity of smooth terms
conCurv = concurvity(redMod,full=FALSE)
round(conCurv$estimate,digits=4)

#             para s(Chl0) s(FSLE0) s(Sal0) s(EKE0) s(SSH0) s(Temp0) s(Slope) s(Aspect)
# para         1  0.0000   0.0000  0.0000  0.0000  0.0000   0.0000   0.0000    0.0000
# s(Chl0)      0  1.0000   0.0322  0.3345  0.0776  0.2111   0.2806   0.1455    0.0477
# s(FSLE0)     0  0.0206   1.0000  0.1710  0.0812  0.1191   0.0448   0.0450    0.0267
# s(Sal0)      0  0.1286   0.1268  1.0000  0.1797  0.2476   0.2139   0.1521    0.1053
# s(EKE0)      0  0.0315   0.1023  0.2900  1.0000  0.1574   0.0864   0.0335    0.0395
# s(SSH0)      0  0.1551   0.1434  0.5754  0.2027  1.0000   0.2993   0.2145    0.1298
# s(Temp0)     0  0.2332   0.0840  0.3741  0.1160  0.2677   1.0000   0.0784    0.0809
# s(Slope)     0  0.0856   0.0484  0.2314  0.0551  0.1906   0.0881   1.0000    0.1961
# s(Aspect)    0  0.0361   0.0194  0.1489  0.0348  0.1523   0.0382   0.1538    1.0000

# Sal0 still quite explained by SSH & Temp0, but proceeding anyway

modCompTable = dredge(redMod,
              beta="none",
              evaluate=TRUE,
              trace=TRUE)

# run optimal model
optMod = get.models(modCompTable,subset=1)

# check p-values
PV = summary(optMod$'256')

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(Aspect, bs = "cc", k = 5) + s(Chl0, bs = "cs", k = 5) + 
#   s(EKE0, bs = "cs", k = 5) + s(FSLE0, bs = "cs", k = 5) + 
#   s(Sal0, bs = "cs", k = 4) + s(Slope, bs = "cc", k = 5) + 
#   s(SSH0, bs = "cs", k = 5) + s(Temp0, bs = "cs", k = 5) + 
#   1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 3.263605   0.002586    1262   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df   Chi.sq p-value    
# s(Aspect) 2.983      3  2966.35  <2e-16 ***
#   s(Chl0)   3.979      4  4500.43  <2e-16 ***
#   s(EKE0)   3.905      4    56.59  <2e-16 ***
#   s(FSLE0)  3.950      4  1910.39  <2e-16 ***
#   s(Sal0)   2.867      3  7528.63  <2e-16 ***
#   s(Slope)  2.994      3  1886.67  <2e-16 ***
#   s(SSH0)   3.986      4 15254.86  <2e-16 ***
#   s(Temp0)  3.999      4 10570.58  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.567   Deviance explained =   60%
# UBRE = 17.128  Scale est. = 1         n = 10477


# plot
plot.gam(optMod$'256',all.terms=TRUE)



