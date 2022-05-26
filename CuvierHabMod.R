library(tidyverse)
# library(car)
library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

## GAM approach ---------------------
# Regional model
spec = 'Cuvier'
outDir = "E:/ModelingCovarData/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/UD28_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

# Check which covars are correlated w presence to determine starting covar list

# Test for how a term should be included in the model
smoothVarList = c ("FSLE0",
                  "Sal0",
                  "Sal700",
                  "SSH",
                  "Temp0",
                  "Temp700",
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

          # linMod             threeKnots         fourKnots          fiveKnots          Best       
# FSLE0   "520088.534596979" "502724.700278472" "502631.878618437" "500301.060225707" "fiveKnots"
# Sal0    "486235.085207736" "411360.941415359" "387124.722844731" "387743.671350693" "fourKnots"
# Sal700  "464609.314309435" "408543.686862355" "387910.578345626" "387757.127764851" "fiveKnots"
# SSH     "315223.54295346"  "313716.500030996" "310937.159438105" "309552.283203949" "fiveKnots"
# Temp0   "376639.997067505" "346409.991686607" "345047.680730214" "344559.370920584" "fiveKnots"
# Temp700 "403429.944778358" "371727.212253311" "359909.516906043" "357272.95604992"  "fiveKnots"
# Slope   "471243.092863463" "378355.828929576" "378355.828929576" "378321.313829823" "fiveKnots"
# Aspect  "519177.311573147" "518337.881472128" "518337.881472128" "512187.826366409" "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=5)
              + s(Slope,bs="cc",k=5)
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
# s(FSLE0)      0   1.0000  0.1710    0.1125  0.1191   0.0448     0.1255   0.0450    0.0267
# s(Sal0)       0   0.1268  1.0000    0.6319  0.2476   0.2139     0.5615   0.1521    0.1053
# s(Sal700)     0   0.1531  0.8001    1.0000  0.3452   0.2256     0.6636   0.1851    0.1746
# s(SSH0)       0   0.1434  0.5754    0.3864  1.0000   0.2993     0.4568   0.2145    0.1298
# s(Temp0)      0   0.0840  0.3741    0.3292  0.2677   1.0000     0.3941   0.0784    0.0809
# s(Temp700)    0   0.1492  0.7265    0.6841  0.3361   0.3469     1.0000   0.1491    0.1256
# s(Slope)      0   0.0484  0.2314    0.1767  0.1906   0.0881     0.1758   1.0000    0.1961
# s(Aspect)     0   0.0194  0.1489    0.1216  0.1523   0.0382     0.0618   0.1538    1.0000


# Sal0 problematic w SSH, Sal700, Temp0, less w Chl0
# Sal700 prolematic w Sal0, Temp700
# Taking out Sal0, Sal700, SSH

# run reduced model
redMod = gam(Pres ~ s(FSLE0,bs="cs",k=5)
             # + s(Sal0,bs="cs",k=4)
             # + s(Sal700,bs="cs",k=5)
             # + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5)
             + s(Temp700,bs="cs",k=5)
             + s(Slope,bs="cc",k=5)
             + s(Aspect,bs="cc",k=5),
             data=data,
             family=poisson,
             gamma=1.4,
             na.action="na.fail")


# check concurvity of smooth terms
conCurv = concurvity(redMod,full=FALSE)
round(conCurv$estimate,digits=4)

#             para s(FSLE0) s(Temp0) s(Temp700) s(Slope) s(Aspect)
# para          1   0.0000   0.0000     0.0000   0.0000    0.0000
# s(FSLE0)      0   1.0000   0.0448     0.1255   0.0450    0.0267
# s(Temp0)      0   0.0840   1.0000     0.3941   0.0784    0.0809
# s(Temp700)    0   0.1492   0.3469     1.0000   0.1491    0.1256
# s(Slope)      0   0.0484   0.0881     0.1758   1.0000    0.1961
# s(Aspect)     0   0.0194   0.0382     0.0618   0.1538    1.0000

# Sal700 still quite explained by SSH, but proceeding anyway

modCompTable = dredge(redMod,
                      beta="none",
                      evaluate=TRUE,
                      trace=TRUE)

# run optimal model
optMod = get.models(modCompTable,subset=1)
optMod = optMod[[names(optMod)]]

# check p-values
PV = summary(optMod)$s.pv
summary(optMod)

# Family: poisson
# Link function: log
# 
# Formula:
#   Pres ~ s(Aspect, bs = "cc", k = 5) + s(FSLE0, bs = "cs",
#                                          k = 5) + s(Slope, bs = "cc", k = 5) + s(Temp0, bs = "cs",
#                                                                                  k = 5) + s(Temp700, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) 3.298419   0.002507    1316   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(Aspect)  2.999      3   4714  <2e-16 ***
#   s(FSLE0)   3.961      4   6894  <2e-16 ***
#   s(Slope)   2.996      3   8245  <2e-16 ***
#   s(Temp0)   3.998      4  46248  <2e-16 ***
#   s(Temp700) 3.995      4  14260  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.528   Deviance explained = 55.1%
# UBRE = 19.376  Scale est. = 1         n = 10477


# plot
png(filename=paste(outDir,'/',spec,'/','Cuvier_allSites.png',sep=""),width=600,height=600)
plot.gam(optMod,all.terms=TRUE,rug=TRUE,scale=0,pages=1)
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
    # same terms as redMod from regional model, but no slope or aspect
    fullSiteMod = gam(Pres ~ s(FSLE0,bs="cs",k=5)
                      # + s(Sal0,bs="cs",k=4)
                      # + s(Sal700,bs="cs",k=5)
                      # + s(SSH0,bs="cs",k=5)
                      + s(Temp0,bs="cs",k=5)
                      + s(Temp700,bs="cs",k=5),
                      # + s(Slope,bs="cc",k=5)
                      # + s(Aspect,bs="cc",k=5),
                      data=siteData,
                      family=poisson,
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
        # update model
        redMod<-eval(parse(text=paste("update(optSiteMod, . ~ . - ", paste(badVars,collapse="-"), ")", sep="")))
        redSitePV = summary(redMod)$s.pv
        if (!any(redSitePV>=0.05)){
          siteModList[[sites[i]]] = redMod
          pValList[[sites[i]]] = redSitePV
          flag=0
        }
      }
    } else {
      siteModList[[sites[i]]] = optSiteMod
      pValList[[sites[i]]] = sitePV
    }
    
    png(filename=paste(outDir,'/',spec,'/','Cuvier_',sites[i],'.png',sep=""),width=600,height=600)
    plot.gam(siteModList[[sites[i]]],all.terms=TRUE,rug=TRUE,scale=0,pages=1)
    while (dev.cur()>1) {dev.off()}
    
  }
}

save(siteModList,pValList,siteModCompList,file=paste(outDir,'/',spec,'/','SiteSpecificModels.Rdata',sep=""))

