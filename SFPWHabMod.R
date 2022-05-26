library(tidyverse)
library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

## GAM approach ---------------------
# Regional model
spec = 'SFPW'
outDir = "E:/ModelingCovarData/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/UD26_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

# Check which covars are correlated w presence to determine starting covar list

# Test for how a term should be included in the model
smoothVarList = c ("Chl0",
                   "FSLE0",
                   "Sal0",
                   "Sal700",
                   "SSH0",
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

#             linMod             threeKnots         fourKnots          fiveKnots          Best       
# Chl0    "173931.539455583" "164401.90983094"  "156155.115172241" "156013.261499281" "fiveKnots"
# FSLE0   "173741.349301875" "170972.349325979" "165872.994538534" "163957.516696382" "fiveKnots"
# Sal0    "177870.753516572" "177199.035056446" "172632.17245835"  "171558.308456471" "fiveKnots"
# Sal700  "178177.295450721" "174602.559816833" "167903.415368144" "162910.972607726" "fiveKnots"
# SSH0    "167877.347367725" "135603.381048807" "135107.411940843" "135184.500097775" "fourKnots"
# Temp0   "177924.455758584" "166214.091499445" "165995.039170708" "165416.245772813" "fiveKnots"
# Temp700 "177740.652544819" "166433.625252021" "161304.991843192" "161409.92895673"  "fourKnots"
# Slope   "163692.761760395" "134021.328642419" "134021.328642419" "127361.218266489" "fiveKnots"
# Aspect  "160392.189591788" "148968.245769183" "148968.245769183" "143476.574508039" "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
              +  s(FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=4)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=4)
              + s(Slope,bs="cc",k=5)
              + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail")

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)

#             para s(Chl0) s(FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(Slope) s(Aspect)
# para          1  0.0000   0.0000  0.0000    0.0000  0.0000   0.0000     0.0000   0.0000    0.0000
# s(Chl0)       0  1.0000   0.0322  0.3345    0.2248  0.6304   0.2806     0.4450   0.1455    0.0477
# s(FSLE0)      0  0.0206   1.0000  0.1710    0.1125  0.1557   0.0448     0.2162   0.0450    0.0267
# s(Sal0)       0  0.1286   0.1268  1.0000    0.6319  0.8188   0.2139     0.8140   0.1521    0.1053
# s(Sal700)     0  0.1230   0.1531  0.8001    1.0000  0.8437   0.2256     0.8758   0.1851    0.1746
# s(SSH0)       0  0.1525   0.1151  0.5624    0.3831  1.0000   0.3035     0.7746   0.1911    0.1320
# s(Temp0)      0  0.2332   0.0840  0.3741    0.3292  0.5223   1.0000     0.5096   0.0784    0.0809
# s(Temp700)    0  0.1227   0.1490  0.7181    0.6400  0.8519   0.3353     1.0000   0.1468    0.0996
# s(Slope)      0  0.0856   0.0484  0.2314    0.1767  0.4729   0.0881     0.3305   1.0000    0.1961
# s(Aspect)     0  0.0361   0.0194  0.1489    0.1216  0.1070   0.0382     0.1043   0.1538    1.0000

# Sal0 problematic w Chl0, Sal700, SSH, Temp700, less with temp0
# Sal700 prolematic w Sal0, Temp700
# SSH problematic with Chl0, Sal0, Sal700, Temp0, Temp700
# Temp700 problematic with Sal0, Sal700, SSH, Temp0
# Taking out SSH, Sal0, Temp700

# run reduced model
redMod =gam(Pres ~ s(Chl0,bs="cs",k=5)
            +  s(FSLE0,bs="cs",k=5)
            # + s(Sal0,bs="cs",k=4)
            + s(Sal700,bs="cs",k=5)
            # + s(SSH0,bs="cs",k=4)
            + s(Temp0,bs="cs",k=5)
            # + s(Temp700,bs="cs",k=4)
            + s(Slope,bs="cc",k=5)
            + s(Aspect,bs="cc",k=5),
             data=data,
             family=poisson,
             gamma=1.4,
             na.action="na.fail")


# check concurvity of smooth terms
conCurv = concurvity(redMod,full=FALSE)
round(conCurv$estimate,digits=4)

#             para s(Chl0) s(FSLE0) s(Sal700) s(Temp0) s(Slope) s(Aspect)
# para         1  0.0000   0.0000    0.0000   0.0000   0.0000    0.0000
# s(Chl0)      0  1.0000   0.0322    0.2248   0.2806   0.1455    0.0477
# s(FSLE0)     0  0.0206   1.0000    0.1125   0.0448   0.0450    0.0267
# s(Sal700)    0  0.1230   0.1531    1.0000   0.2256   0.1851    0.1746
# s(Temp0)     0  0.2332   0.0840    0.3292   1.0000   0.0784    0.0809
# s(Slope)     0  0.0856   0.0484    0.1767   0.0881   1.0000    0.1961
# s(Aspect)    0  0.0361   0.0194    0.1216   0.0382   0.1538    1.0000

# Sal700 still quite explained by Temp0, but proceeding anyway

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
#   Pres ~ s(Aspect, bs = "cc", k = 5) + s(Chl0, bs = "cs", 
#                                          k = 5) + s(FSLE0, bs = "cs", k = 5) + s(Sal700, bs = "cs", 
#                                                                                  k = 5) + s(Slope, bs = "cc", k = 5) + s(Temp0, bs = "cs", 
#                                                                                                                          k = 5) + 1
# 
# Parametric coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 0.817079   0.009114   89.65   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(Aspect) 2.997      3 6071.1  <2e-16 ***
#   s(Chl0)   4.000      4 3153.7  <2e-16 ***
#   s(FSLE0)  3.982      4  513.2  <2e-16 ***
#   s(Sal700) 3.919      4  703.7  <2e-16 ***
#   s(Slope)  2.994      3 5973.5  <2e-16 ***
#   s(Temp0)  3.986      4 8303.5  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.415   Deviance explained = 49.6%
# UBRE = 6.7549  Scale est. = 1         n = 10477


# plot
png(filename=paste(outDir,'/',spec,'/','SFPW_allSites.png',sep=""),width=600,height=600)
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
    fullSiteMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
                      +  s(FSLE0,bs="cs",k=5)
                      # + s(Sal0,bs="cs",k=4)
                      + s(Sal700,bs="cs",k=5)
                      # + s(SSH0,bs="cs",k=4)
                      + s(Temp0,bs="cs",k=5),
                      # + s(Temp700,bs="cs",k=4)
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
    
    png(filename=paste(outDir,'/',spec,'/','SFPW_',sites[i],'.png',sep=""),width=600,height=600)
    plot.gam(siteModList[[sites[i]]],all.terms=TRUE,rug=TRUE,scale=0,pages=1)
    while (dev.cur()>1) {dev.off()}
    
  }
}

save(siteModList,pValList,siteModCompList,file=paste(outDir,'/',spec,'/','SiteSpecificModels.Rdata',sep=""))

