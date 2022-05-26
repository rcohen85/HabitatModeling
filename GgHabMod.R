library(tidyverse)
library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

# library(tidyverse)
# library(splines2)
# library(geepack)
# source("getPvalues.R")

data = data.frame(read.csv('J:/Chpt_3/ModelData/Risso_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

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

## GAM approach ---------------------
# Regional model
spec = 'Risso'
outDir = "J:/Chpt_3/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("Chl0",
                  "FSLE0",
                  "Sal400",
                  "SSH0",
                  "Temp0",
                  "VelMag400",
                  "Slope",
                  "Aspect")

# Test for how a term should be included in the model
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
# Chl0      "261312.471668402" "252419.905096666" "240554.844722747" "237682.086219538" "fiveKnots"
# FSLE0     "264553.795241769" "260440.831334124" "255574.09133929"  "255267.899784281" "fiveKnots"
# Sal400    "241644.107670266" "198934.166020971" "192674.656880729" "192320.527800769" "fiveKnots"
# SSH0      "196025.04362009"  "195910.540107404" "195895.557766024" "194503.345499942" "fiveKnots"
# Temp0     "252833.460857668" "236419.316968059" "235034.561409143" "234993.26356563"  "fiveKnots"
# VelMag400 "253208.581673729" "251608.011768822" "251356.27842725"  "250762.422753389" "fiveKnots"
# Slope     "273274.371807135" "250070.644800729" "236963.97047091"  "229526.092130719" "fiveKnots"
# Aspect    "220977.026245278" "261010.741102481" "261010.741102481" "255274.826466453" "linMod"  

# run full model
fullMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
              + s(FSLE0,bs="cs",k=5)
              + s(Sal400,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(VelMag400,bs="cc",k=5)
              + s(Slope,bs="cs",k=5)
              + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
              method="REML",
              select=TRUE,
              gamma=1.4,
              na.action="na.fail")

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)

#               para s(Chl0) s(FSLE0) s(Sal400) s(SSH0) s(Temp0) s(VelMag400) s(Slope) s(Aspect)
# para            1  0.0000   0.0000    0.0000  0.0000   0.0000       0.0000   0.0000    0.0000
# s(Chl0)         0  1.0000   0.0322    0.2505  0.2111   0.2806       0.0397   0.1385    0.0477
# s(FSLE0)        0  0.0206   1.0000    0.1212  0.1191   0.0448       0.1087   0.0401    0.0267
# s(Sal400)       0  0.1225   0.1590    1.0000  0.3233   0.2339       0.2634   0.1887    0.1800
# s(SSH0)         0  0.1551   0.1434    0.4346  1.0000   0.2993       0.2965   0.2119    0.1298
# s(Temp0)        0  0.2332   0.0840    0.3260  0.2677   1.0000       0.1271   0.0717    0.0809
# s(VelMag400)    0  0.0197   0.1130    0.1788  0.2364   0.0806       1.0000   0.0194    0.1032
# s(Slope)        0  0.0919   0.0577    0.2062  0.2844   0.0945       0.0786   1.0000    0.2540
# s(Aspect)       0  0.0361   0.0194    0.1240  0.1523   0.0382       0.1574   0.1911    1.0000

# not removing anything

modCompTable = dredge(fullMod,
                      beta="none",
                      evaluate=TRUE,
                      trace=TRUE)

# run optimal model
optMod = get.models(modCompTable,subset=1)
optMod = optMod[[names(optMod)]]
save(optMod,modCompTable,file=paste(outDir,'/',spec,'/','RegionalModel.Rdata',sep=""))

# check p-values
PV = summary(optMod)$s.pv
summary(optMod)

# Family: poisson 
# Link function: log 
# 
# Formula:
#   Pres ~ s(Aspect, bs = "cc", k = 5) + s(Chl0, bs = "cs", k = 5) + 
#   s(FSLE0, bs = "cs", k = 5) + s(Sal400, bs = "cs", k = 5) + 
#   s(Slope, bs = "cs", k = 5) + s(SSH0, bs = "cs", k = 5) + 
#   s(Temp0, bs = "cs", k = 5) + s(VelMag400, bs = "cc", k = 5) + 
#   1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 1.303678   0.009342   139.5   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq p-value    
# s(Aspect)    2.999      3 9377.78  <2e-16 ***
#   s(Chl0)      3.817      4 1427.57  <2e-16 ***
#   s(FSLE0)     3.930      4 1630.82  <2e-16 ***
#   s(Sal400)    3.993      4 2287.23  <2e-16 ***
#   s(Slope)     3.963      4 1393.48  <2e-16 ***
#   s(SSH0)      3.992      4 8885.22  <2e-16 ***
#   s(Temp0)     3.960      4 3268.27  <2e-16 ***
#   s(VelMag400) 2.371      3   37.75  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.302   Deviance explained = 47.2%
# -REML =  55987  Scale est. = 1         n = 10477


# plot
png(filename=paste(outDir,'/',spec,'/',spec,'_allSites.png',sep=""),width=600,height=600)
plot.gam(optMod,all.terms=TRUE,rug=TRUE,pages=1,scale=0)
while (dev.cur()>1) {dev.off()}


# Site-specific models ---------------
sites = unique(data$Site)
siteModList = list()
pValList = list()
siteModCompList = list()
for (i in 1:length(sites)){
  
  siteInd = which(!is.na(str_match(data$Site,sites[i])))
  siteData = data[siteInd,]
  
  if (sum(siteData$Pres>0)>25){
    
    fullSiteMod = gam(Pres ~ s(Chl0,bs="cs",k=5) # same terms as redMod from regional model, but no slope or aspect
                      + s(FSLE0,bs="cs",k=5)
                      + s(Sal400,bs="cs",k=5)
                      + s(SSH0,bs="cs",k=5)
                      + s(Temp0,bs="cs",k=5)
                      + s(VelMag400,bs="cc",k=5),
                      data=siteData,
                      family=poisson,
                      method="REML",
                      select=TRUE,
                      gamma=1.4,
                      na.action="na.fail")
    
    siteModCompTable = dredge(fullSiteMod,
                              beta="none",
                              evaluate=TRUE,
                              trace=TRUE)
    siteModCompList[[sites[i]]] = siteModCompTable
    
    optSiteMod = get.models(siteModCompTable,subset=1)
    optSiteMod = optSiteMod[[names(optSiteMod)]]
    sitePV = summary(optSiteMod)$s.pv
    
    
    if (any(sitePV>=0.05)){ # Remove non-significant terms & re-run model iteratively until only signif covars remain
      flag = 1
      while (flag==1){
        # get terms from formula as strings
        thisForm = as.character(optSiteMod$formula)[3]
        startSmooth = str_locate_all(thisForm,'s\\(')[[1]][,1]
        termInd = str_locate_all(thisForm,'\\+')[[1]][,1]
        termInd = c(0,termInd,str_length(thisForm)+1)
        allTerms = character()
        for (j in 1:length(termInd)-1){
          thisTerm = str_sub(thisForm,start=termInd[j]+1,end=termInd[j+1]-1)
          allTerms = c(allTerms,thisTerm)
        }
        # identify which terms were non-significant
        badVars = allTerms[sitePV>=0.05]
        dontNeed = which(!is.na(str_match(badVars,"1")))
        badVars = badVars[-dontNeed]
        # update model
        optSiteMod<-eval(parse(text=paste("update(optSiteMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
        sitePV = summary(optSiteMod)$s.pv
        if (!any(sitePV>=0.05)){
          siteModList[[sites[i]]] = optSiteMod
          pValList[[sites[i]]] = sitePV
          flag=0
        }
      }
    } else {
      siteModList[[sites[i]]] = optSiteMod
      pValList[[sites[i]]] = sitePV
    }
    
    png(filename=paste(outDir,'/',spec,'/',spec,'_',sites[i],'.png',sep=""),width=600,height=600)
    plot.gam(siteModList[[sites[i]]],all.terms=TRUE,rug=TRUE,pages=1,main=sites[i],scale=0)
    while (dev.cur()>1) {dev.off()}
    
  }
}

save(siteModList,pValList,siteModCompList,file=paste(outDir,'/',spec,'/','SiteSpecificModels.Rdata',sep=""))

