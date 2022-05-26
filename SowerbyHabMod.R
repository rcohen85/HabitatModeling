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

data = data.frame(read.csv('J:/Chpt_3/ModelData/Sowerby_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

## GAM approach ---------------------
# Regional model
spec = 'Sowerby'
outDir = "J:/Chpt_3/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("Chl0",
                  "FSLE0",
                  "Sal0",
                  "Sal700",
                  "SSH0",
                  "Temp0",
                  "Temp700",
                  "VelMag0",
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
# Chl0    "35285.8070203761" "33496.1306791603" "31919.3472884446" "31887.696428841"  "fiveKnots"
# FSLE0   "36058.0740385027" "36022.7297613808" "35824.2365576127" "35807.5206416723" "fiveKnots"
# Sal0    "32876.9400761289" "30727.147217092"  "29402.5332038039" "29408.2490413506" "fourKnots"
# Sal700  "32849.3418575005" "31084.3627676835" "29229.9143121705" "29182.4233706146" "fiveKnots"
# SSH0    "29449.7147840002" "29052.9833044351" "29047.9460718174" "29029.0085640016" "fiveKnots"
# Temp0   "33027.7144177467" "32597.1114003647" "32507.1552647256" "32447.6185267845" "fiveKnots"
# Temp700 "30347.1300806573" "29347.1594844628" "29113.377268367"  "29108.0814772288" "fiveKnots"
# VelMag0 "35312.2264463077" "35008.9366624951" "35000.1201284772" "34882.5751749246" "fiveKnots"
# Slope   "33206.7189659358" "29479.3463359707" "29083.0627070293" "28197.9595158366" "fiveKnots"
# Aspect  "36988.897533702"  "34209.7310488775" "34209.7310488775" "34173.639065464"  "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
              + s(FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=5)
              + s(VelMag0,bs="cc",k=5)
              + s(Slope,bs="cs",k=5)
              + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail")

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)

#             para s(Chl0) s(FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(VelMag0) s(Slope) s(Aspect)
# para          1  0.0000   0.0000  0.0000    0.0000  0.0000   0.0000     0.0000     0.0000   0.0000    0.0000
# s(Chl0)       0  1.0000   0.0322  0.3345    0.2248  0.2111   0.2806     0.2390     0.0321   0.1385    0.0477
# s(FSLE0)      0  0.0206   1.0000  0.1710    0.1125  0.1191   0.0448     0.1255     0.1011   0.0401    0.0267
# s(Sal0)       0  0.1286   0.1268  1.0000    0.6319  0.2476   0.2139     0.5615     0.1556   0.1414    0.1053
# s(Sal700)     0  0.1230   0.1531  0.8001    1.0000  0.3452   0.2256     0.6636     0.2435   0.1839    0.1746
# s(SSH0)       0  0.1551   0.1434  0.5754    0.3864  1.0000   0.2993     0.4568     0.2710   0.2119    0.1298
# s(Temp0)      0  0.2332   0.0840  0.3741    0.3292  0.2677   1.0000     0.3941     0.1157   0.0717    0.0809
# s(Temp700)    0  0.1252   0.1492  0.7265    0.6841  0.3361   0.3469     1.0000     0.1731   0.1361    0.1256
# s(VelMag0)    0  0.0173   0.1065  0.2283    0.1490  0.2303   0.0731     0.1344     1.0000   0.0194    0.0932
# s(Slope)      0  0.0919   0.0577  0.2556    0.1940  0.2844   0.0945     0.1849     0.0707   1.0000    0.2540
# s(Aspect)     0  0.0361   0.0194  0.1489    0.1216  0.1523   0.0382     0.0618     0.1415   0.1911    1.0000

# Sal0 problematic w Sal700, SSH, Temp700, less so w Temp0
# Sal700 problematic w Sal0, Temp700, less so w SSH, Temp0
# removing Sal0, Temp700

redMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
              + s(FSLE0,bs="cs",k=5)
              # + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              # + s(Temp700,bs="cs",k=5)
              + s(VelMag0,bs="cc",k=5)
              + s(Slope,bs="cs",k=5)
              + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
             method="REML",
             select=TRUE,
              gamma=1.4,
              na.action="na.fail")

# check concurvity of smooth terms
conCurv = concurvity(redMod,full=FALSE)
round(conCurv$estimate,digits=4)

#             para s(Chl0) s(FSLE0) s(Sal700) s(SSH0) s(Temp0) s(VelMag0) s(Slope) s(Aspect)
# para          1  0.0000   0.0000    0.0000  0.0000   0.0000     0.0000   0.0000    0.0000
# s(Chl0)       0  1.0000   0.0322    0.2248  0.2111   0.2806     0.0321   0.1385    0.0477
# s(FSLE0)      0  0.0206   1.0000    0.1125  0.1191   0.0448     0.1011   0.0401    0.0267
# s(Sal700)     0  0.1230   0.1531    1.0000  0.3452   0.2256     0.2435   0.1839    0.1746
# s(SSH0)       0  0.1551   0.1434    0.3864  1.0000   0.2993     0.2710   0.2119    0.1298
# s(Temp0)      0  0.2332   0.0840    0.3292  0.2677   1.0000     0.1157   0.0717    0.0809
# s(VelMag0)    0  0.0173   0.1065    0.1490  0.2303   0.0731     1.0000   0.0194    0.0932
# s(Slope)      0  0.0919   0.0577    0.1940  0.2844   0.0945     0.0707   1.0000    0.2540
# s(Aspect)     0  0.0361   0.0194    0.1216  0.1523   0.0382     0.1415   0.1911    1.0000

modCompTable = dredge(redMod,
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
#   s(FSLE0, bs = "cs", k = 5) + s(Sal700, bs = "cs", k = 5) + 
#   s(Slope, bs = "cs", k = 5) + s(SSH0, bs = "cs", k = 5) + 
#   s(Temp0, bs = "cs", k = 5) + s(VelMag0, bs = "cc", k = 5) + 
#   1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -6.5055     0.6377   -10.2   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq p-value    
# s(Aspect)  2.976      3  769.66  <2e-16 ***
#   s(Chl0)    2.300      4  120.99  <2e-16 ***
#   s(FSLE0)   2.878      4   36.98  <2e-16 ***
#   s(Sal700)  3.968      4  240.11  <2e-16 ***
#   s(Slope)   3.938      4 3255.80  <2e-16 ***
#   s(SSH0)    3.097      4  356.95  <2e-16 ***
#   s(Temp0)   3.852      4  122.61  <2e-16 ***
#   s(VelMag0) 2.445      3   51.77  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.315   Deviance explained = 50.7%
# -REML = 7925.6  Scale est. = 1         n = 10477


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
                          + s(Sal700,bs="cs",k=5)
                          + s(SSH0,bs="cs",k=5)
                          + s(Temp0,bs="cs",k=5)
                          + s(VelMag0,bs="cc",k=5),
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

