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

data = data.frame(read.csv('J:/Chpt_3/ModelData/True_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

## GAM approach ---------------------
# Regional model
spec = 'True'
outDir = "J:/Chpt_3/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

# Check which covars are correlated w presence to determine starting covar list
smoothVarList = c("FSLE0",
                  "Sal0",
                  "Sal700",
                  "SSH0",
                  "Temp0",
                  "Temp700",
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

#           linMod             threeKnots         fourKnots          fiveKnots          Best       
# FSLE0   "21064.6398221111" "21063.7550941797" "21065.7892192821" "21054.8712017735" "fiveKnots"
# Sal0    "20591.8886962384" "19508.3247970513" "19036.7101651424" "19025.2679925543" "fiveKnots"
# Sal700  "20562.1233107988" "19144.1219894864" "18796.3368344045" "18785.715569249"  "fiveKnots"
# SSH0    "18652.8959407941" "18550.9983715304" "18335.7194326727" "18303.8783518756" "fiveKnots"
# Temp0   "20310.3205801655" "19807.1689108074" "19778.1605284679" "19707.5989740754" "fiveKnots"
# Temp700 "19558.915555238"  "18543.3912988503" "18546.5624987243" "18454.6488797136" "fiveKnots"
# Slope   "21982.4332449449" "20162.3313689192" "19959.4713055543" "19959.0459960718" "fiveKnots"
# Aspect  "20794.3295873248" "22179.3729901377" "22179.3729901377" "22145.1072714889" "linMod"

# run full model
fullMod = gam(Pres ~ s(FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=5)
              + s(Slope,bs="cs",k=5)
              + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail")

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)

#             para s(FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(Slope) s(Aspect)
# para          1   0.0000  0.0000    0.0000  0.0000   0.0000     0.0000   0.0000    0.0000
# s(FSLE0)      0   1.0000  0.1162    0.1125  0.1191   0.0448     0.1255   0.0401    0.0267
# s(Sal0)       0   0.1560  1.0000    0.6640  0.2913   0.2460     0.5787   0.1452    0.1161
# s(Sal700)     0   0.1531  0.6765    1.0000  0.3452   0.2256     0.6636   0.1839    0.1746
# s(SSH0)       0   0.1434  0.4231    0.3864  1.0000   0.2993     0.4568   0.2119    0.1298
# s(Temp0)      0   0.0840  0.3360    0.3292  0.2677   1.0000     0.3941   0.0717    0.0809
# s(Temp700)    0   0.1492  0.5950    0.6841  0.3361   0.3469     1.0000   0.1361    0.1256
# s(Slope)      0   0.0577  0.1944    0.1940  0.2844   0.0945     0.1849   1.0000    0.2540
# s(Aspect)     0   0.0194  0.1047    0.1216  0.1523   0.0382     0.0618   0.1911    1.0000

# Sal0 problematic w Sal700, SSH, Temp700
# Sal700 problematic w Sal0, Temp700
# Taking out Sal0 and Temp700

redMod = gam(Pres ~ s(FSLE0,bs="cs",k=5)
             # + s(Sal0,bs="cs",k=5)
             + s(Sal700,bs="cs",k=5)
             + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5)
             # + s(Temp700,bs="cs",k=5)
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

#             para s(FSLE0) s(Sal700) s(SSH0) s(Temp0) s(Slope) s(Aspect)
# para         1   0.0000    0.0000  0.0000   0.0000   0.0000    0.0000
# s(FSLE0)     0   1.0000    0.1125  0.1191   0.0448   0.0401    0.0267
# s(Sal700)    0   0.1531    1.0000  0.3452   0.2256   0.1839    0.1746
# s(SSH0)      0   0.1434    0.3864  1.0000   0.2993   0.2119    0.1298
# s(Temp0)     0   0.0840    0.3292  0.2677   1.0000   0.0717    0.0809
# s(Slope)     0   0.0577    0.1940  0.2844   0.0945   1.0000    0.2540
# s(Aspect)    0   0.0194    0.1216  0.1523   0.0382   0.1911    1.0000

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
#   Pres ~ s(Aspect, bs = "cc", k = 5) + s(FSLE0, bs = "cs", k = 5) + 
#   s(Sal700, bs = "cs", k = 5) + s(Slope, bs = "cs", k = 5) + 
#   s(SSH0, bs = "cs", k = 5) + s(Temp0, bs = "cs", k = 5) + 
#   1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.47711    0.09224  -26.86   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(Aspect) 2.960      3 409.28  <2e-16 ***
#   s(FSLE0)  3.795      4 220.88  <2e-16 ***
#   s(Sal700) 1.733      4  44.64  <2e-16 ***
#   s(Slope)  3.832      4 144.13  <2e-16 ***
#   s(SSH0)   3.788      4 335.69  <2e-16 ***
#   s(Temp0)  3.647      4 169.79  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0981   Deviance explained = 26.9%
# -REML =   6209  Scale est. = 1         n = 10477


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
    
    fullSiteMod = gam(Pres ~ s(FSLE0,bs="cs",k=5) # same terms as redMod from regional model, but no slope or aspect
                      + s(Sal700,bs="cs",k=5)
                      + s(SSH0,bs="cs",k=5)
                      + s(Temp0,bs="cs",k=5),
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

