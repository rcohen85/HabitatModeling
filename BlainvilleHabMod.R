library(tidyverse)
# library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

## GAM approach ---------------------
# Regional model
spec = 'Blainville'
outDir = "E:/ModelingCovarData/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/Blainville_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

# Check which covars are correlated w presence to determine starting covar list

# Test for how a term should be included in the model
smoothVarList = c ("EKE0",
                   "Sal0",
                   "Sal700",
                   "Temp0",
                   "VelAsp700",
                   "VelMag700",
                   "Slope",
                   "Aspect")


modOpts = c("linMod","threeKnots","fourKnots","fiveKnots")
AIC_votes = matrix(nrow=length(smoothVarList),ncol=5)

for (i in 1:(length(smoothVarList))){
  
  if (str_detect(smoothVarList[i],"Asp")){
    bs = "cc"
  } else { bs = "cs"}
  
  modelCall = paste('gam(Pres~data$',smoothVarList[i],',data=data,family=poisson,method="REML",select=TRUE)',sep="")
  linMod = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=3),data=data,family=poisson,method="REML",select=TRUE)',sep="")
  smoothMod1 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=4),data=data,family=poisson,method="REML",select=TRUE)',sep="")
  smoothMod2 = eval(parse(text=modelCall))
  
  modelCall = paste('gam(Pres~s(data$',smoothVarList[i],',bs="',bs,'",k=5),data=data,family=poisson,method="REML",select=TRUE)',sep="")
  smoothMod3 = eval(parse(text=modelCall))
  
  AIC_votes[i,1:4] = c(AIC(linMod)[[1]],AIC(smoothMod1)[[1]],AIC(smoothMod2)[[1]],AIC(smoothMod3)[[1]])
  AIC_votes[i,5] = modOpts[which.min(AIC_votes[i,1:4])]
}

colnames(AIC_votes) = c(modOpts,"Best")
rownames(AIC_votes) = smoothVarList[]
AIC_votes

#               linMod             threeKnots         fourKnots          fiveKnots          Best        
# EKE0      "50390.2512384755" "49273.5908057239" "49016.5298778083" "47849.1000439767" "fiveKnots" 
# Sal0      "27042.7451409243" "26963.1097836539" "26964.0463900206" "26966.9033951201" "threeKnots"
# Sal700    "25099.894702228"  "21840.011549816"  "21839.9795669367" "21839.9256211684" "fiveKnots" 
# Temp0     "44910.8032564516" "41679.9956984451" "40693.3895658277" "40417.0659111374" "fiveKnots" 
# VelAsp700 "50110.0867110415" "49585.0981321812" "49585.0981321812" "45286.7430131016" "fiveKnots" 
# VelMag700 "48944.1798858157" "47942.7067631457" "47842.4220947324" "47819.3618459435" "fiveKnots" 
# Slope     "36549.8039787732" "13226.6675004817" "13127.2732566947" "13126.1841830934" "fiveKnots" 
# Aspect    "14628.0060309852" "25237.7576542754" "25237.7576542754" "24593.19222504"   "linMod"

# run full model
fullMod = gam(Pres ~ s(EKE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=3)
              + s(Sal700,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(VelAsp700,bs="cc",k=5)
              + s(VelMag700,bs="cs",k=5)
              + s(Slope,bs="cs",k=5)
              + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              select=TRUE)

# check concurvity of smooth terms
conCurv = concurvity(fullMod,full=FALSE)
round(conCurv$estimate,digits=4)

#               para s(EKE0) s(Sal0) s(Sal700) s(Temp0) s(VelAsp700) s(VelMag700) s(Slope) s(Aspect)
# para            1  0.0000  0.0000    0.0000   0.0000       0.0000       0.0000   0.0000    0.0000
# s(EKE0)         0  1.0000  0.3109    0.1849   0.0864       0.4787       0.2850   0.0278    0.0395
# s(Sal0)         0  0.1767  1.0000    0.5828   0.2018       0.1750       0.1940   0.1184    0.0794
# s(Sal700)       0  0.2105  0.8715    1.0000   0.2256       0.2080       0.3313   0.1839    0.1746
# s(Temp0)        0  0.1160  0.4218    0.3292   1.0000       0.1069       0.1582   0.0717    0.0809
# s(VelAsp700)    0  0.4829  0.3620    0.2078   0.0901       1.0000       0.2578   0.0344    0.0414
# s(VelMag700)    0  0.1767  0.2672    0.1901   0.0891       0.1503       1.0000   0.0221    0.1093
# s(Slope)        0  0.0649  0.3808    0.1940   0.0945       0.0788       0.0759   1.0000    0.2540
# s(Aspect)       0  0.0348  0.1343    0.1216   0.0382       0.0493       0.1946   0.1911    1.0000

# EKE0 problematic with VelAsp700
# Sal0 problematic w Sal700, less w Temp0
# Sal700 prolematic w Sal0, less w Temp0
# Taking out Sal0, EKE0

# run reduced model
redMod = gam(Pres ~ 
               # s(EKE0,bs="cs",k=5)
              # + s(Sal0,bs="cs",k=3)
               s(Sal700,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(VelAsp700,bs="cc",k=5)
              # + s(VelMag700,bs="cs",k=5)
              + s(Slope,bs="cs",k=5)
              + s(Aspect,bs="cc",k=5),
              data=data,
              family=poisson,
              gamma=1.4,
              na.action="na.fail",
              method="REML",
              select=TRUE)

# check concurvity of smooth terms
conCurv = concurvity(redMod,full=FALSE)
round(conCurv$estimate,digits=4)

#               para s(Sal700) s(Temp0) s(VelAsp700) s(Slope) s(Aspect)
# para            1    0.0000   0.0000       0.0000   0.0000    0.0000
# s(Sal700)       0    1.0000   0.2256       0.2080   0.1839    0.1746
# s(Temp0)        0    0.3292   1.0000       0.1069   0.0717    0.0809
# s(VelAsp700)    0    0.2078   0.0901       1.0000   0.0344    0.0414
# s(Slope)        0    0.1940   0.0945       0.0788   1.0000    0.2540
# s(Aspect)       0    0.1216   0.0382       0.0493   0.1911    1.0000

# Sal700 still quite explained by Temp0, but proceeding anyway
# Removed VelMag700 after initial run through bc it was not found to be significant

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
#   Pres ~ s(Aspect, bs = "cc", k = 5) + s(Sal700, bs = "cs", 
#                                          k = 5) + s(Temp0, bs = "cs", k = 5) + s(VelAsp700, 
#                                                                                  bs = "cc", k = 5) + 1
# 
# Parametric coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3362.5      710.7  -4.731 2.23e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                    edf Ref.df  Chi.sq p-value    
#   s(Aspect)    2.921      3 2283.91  <2e-16 ***
#   s(Sal700)    2.923      4   79.98  <2e-16 ***
#   s(Temp0)     2.921      4  441.70  <2e-16 ***
#   s(VelAsp700) 2.949      3  113.79  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.527   Deviance explained = 82.5%
# -REML = 4271.1  Scale est. = 1         n = 10477



# plot
png(filename=paste(outDir,'/',spec,'/','Blainville_allSites.png',sep=""),width=600,height=600)
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
    fullSiteMod = gam(Pres ~ 
                        # s(EKE0,bs="cs",k=5)
                        # + s(Sal0,bs="cs",k=3)
                        s(Sal700,bs="cs",k=5)
                      + s(Temp0,bs="cs",k=5)
                      + s(VelAsp700,bs="cc",k=5),
                      # + s(VelMag700,bs="cs",k=5)
                      # + s(Slope,bs="cs",k=5)
                      # + s(Aspect,bs="cc",k=5),
                      data=siteData,
                      family=poisson,
                      gamma=1.4,
                      na.action="na.fail",
                      method="REML",
                      select=TRUE)
    
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
        if (!is_empty(dontNeed)){
          badVars = badVars[-dontNeed]}
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
    
    png(filename=paste(outDir,'/',spec,'/','Blainville_',sites[i],'.png',sep=""),width=600,height=600)
    plot.gam(siteModList[[sites[i]]],all.terms=TRUE,rug=TRUE,scale=0,pages=1)
    while (dev.cur()>1) {dev.off()}
    
  }
}

save(siteModList,pValList,siteModCompList,file=paste(outDir,'/',spec,'/','SiteSpecificModels.Rdata',sep=""))

