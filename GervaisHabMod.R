library(tidyverse)
# library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

## GAM approach ---------------------
# Regional model
spec = 'Gervais'
outDir = "E:/ModelingCovarData/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/Gervais_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

# Check which covars are correlated w presence to determine starting covar list

# Test for how a term should be included in the model
smoothVarList = c ("Chl0",
                   "FSLE0",
                   "Sal0",
                   "Sal700",
                   "SSH0",
                   "Temp700",
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

#           linMod             threeKnots         fourKnots          fiveKnots          Best        
# Chl0    "118229.927645981" "115489.215748967" "108018.690331556" "106068.782270985" "fiveKnots" 
# FSLE0   "162553.434443166" "155336.305505699" "155286.229403025" "154649.122020783" "fiveKnots" 
# Sal0    "94531.3615646504" "81711.8148979425" "80990.9966251957" "80393.4151280875" "fiveKnots" 
# Sal700  "142256.753203003" "109217.197203478" "107775.905232529" "97193.9785944924" "fiveKnots" 
# SSH0    "109854.618143319" "73579.2486224001" "73352.3604758128" "72859.6904024093" "fiveKnots" 
# Temp700 "134144.502845382" "79052.637656501"  "73326.6538113701" "71534.5753083113" "fiveKnots" 
# Slope   "147405.450630565" "145364.839384447" "134260.094145996" "127176.190477336" "fiveKnots" 
# Aspect  "159900.695276914" "132615.392747505" "132615.392747505" "132637.20029411"  "threeKnots"

# run full model
fullMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
              + s(FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=5)
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

#             para s(Chl0) s(FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp700) s(Slope) s(Aspect)
# para          1  0.0000   0.0000  0.0000    0.0000  0.0000     0.0000   0.0000    0.0000
# s(Chl0)       0  1.0000   0.0322  0.3345    0.1604  0.2111     0.1842   0.1385    0.0477
# s(FSLE0)      0  0.0206   1.0000  0.1710    0.0205  0.1191     0.0321   0.0401    0.0267
# s(Sal0)       0  0.1286   0.1268  1.0000    0.1783  0.2476     0.1928   0.1414    0.1053
# s(Sal700)     0  0.1192   0.0479  0.4151    1.0000  0.2931     0.6843   0.1857    0.0873
# s(SSH0)       0  0.1551   0.1434  0.5754    0.3578  1.0000     0.5074   0.2119    0.1298
# s(Temp700)    0  0.1318   0.0833  0.5454    0.7176  0.4764     1.0000   0.1988    0.1550
# s(Slope)      0  0.0919   0.0577  0.2556    0.3340  0.2844     0.4322   1.0000    0.2540
# s(Aspect)     0  0.0361   0.0194  0.1489    0.1707  0.1523     0.2355   0.1911    1.0000


# Sal0 problematic w SSH, Temp700, less w Chl0, Sal700
# Sal700 problematic w Temp700
# Taking out Temp700, SSH, Sal0

# run reduced model
redMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
             + s(FSLE0,bs="cs",k=5)
             # + s(Sal0,bs="cs",k=4)
             + s(Sal700,bs="cs",k=5)
             # + s(SSH0,bs="cs",k=5)
             # + s(Temp700,bs="cs",k=5)
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

#             para s(Chl0) s(FSLE0) s(Sal700) s(Slope) s(Aspect)
# para         1  0.0000   0.0000    0.0000   0.0000    0.0000
# s(Chl0)      0  1.0000   0.0322    0.1604   0.1385    0.0477
# s(FSLE0)     0  0.0206   1.0000    0.0205   0.0401    0.0267
# s(Sal700)    0  0.1192   0.0479    1.0000   0.1857    0.0873
# s(Slope)     0  0.0919   0.0577    0.3340   1.0000    0.2540
# s(Aspect)    0  0.0361   0.0194    0.1707   0.1911    1.0000

# all looks good!

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
#                                                                                  k = 5) + s(Slope, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -12.641      2.468  -5.123 3.01e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                edf Ref.df Chi.sq p-value    
#   s(Aspect) 2.925      3  219.3  <2e-16 ***
#   s(Chl0)   3.939      4  843.2  <2e-16 ***
#   s(FSLE0)  3.863      4  225.4  <2e-16 ***
#   s(Sal700) 3.184      4  424.9  <2e-16 ***
#   s(Slope)  3.968      4  394.5  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.482   Deviance explained = 73.3%
# -REML =  18849  Scale est. = 1         n = 10477


# plot
png(filename=paste(outDir,'/',spec,'/','Gervais_allSites.png',sep=""),width=600,height=600)
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
                      + s(FSLE0,bs="cs",k=5)
                      # + s(Sal0,bs="cs",k=4)
                      + s(Sal700,bs="cs",k=5),
                      # + s(SSH0,bs="cs",k=5)
                      # + s(Temp700,bs="cs",k=5)
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
    
    png(filename=paste(outDir,'/',spec,'/','Gervais_',sites[i],'.png',sep=""),width=600,height=600)
    plot.gam(siteModList[[sites[i]]],all.terms=TRUE,rug=TRUE,scale=0,pages=1)
    while (dev.cur()>1) {dev.off()}
    
  }
}

save(siteModList,pValList,siteModCompList,file=paste(outDir,'/',spec,'/','SiteSpecificModels.Rdata',sep=""))

