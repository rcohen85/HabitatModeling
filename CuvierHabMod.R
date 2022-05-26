library(tidyverse)
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

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/Cuvier_masterDF.csv'))
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
# FSLE0   "143518.029434583" "139700.990717399" "136581.634889274" "136493.20158728"  "fiveKnots"
# Sal0    "153118.970496731" "153007.155354861" "138251.547146199" "131995.945355285" "fiveKnots"
# Sal700  "153926.405063222" "153828.301284669" "135649.612216234" "129176.277855976" "fiveKnots"
# SSH     "166827.993364552" "95270.4496259113" "83489.9374853719" "82791.4859578095" "fiveKnots"
# Temp0   "158889.846656827" "158876.379796124" "157468.377276248" "157446.412951153" "fiveKnots"
# Temp700 "155067.737181327" "149396.209582999" "137529.18534134"  "138160.689819933" "fourKnots"
# Slope   "119444.909524396" "109800.390547004" "108806.710522833" "103469.192587452" "fiveKnots"
# Aspect  "141184.85109818"  "140609.895994276" "140609.895994276" "134940.422019047" "fiveKnots"

# run full model
fullMod = gam(Pres ~ s(FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=5)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=5)
              + s(Temp700,bs="cs",k=4)
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

#             para s(FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(Slope) s(Aspect)
# para          1   0.0000  0.0000    0.0000  0.0000   0.0000     0.0000   0.0000    0.0000
# s(FSLE0)      0   1.0000  0.1162    0.1125  0.1191   0.0448     0.2162   0.0401    0.0267
# s(Sal0)       0   0.1560  1.0000    0.6640  0.2913   0.2460     0.8290   0.1452    0.1161
# s(Sal700)     0   0.1531  0.6765    1.0000  0.3452   0.2256     0.8758   0.1839    0.1746
# s(SSH0)       0   0.1434  0.4231    0.3864  1.0000   0.2993     0.7826   0.2119    0.1298
# s(Temp0)      0   0.0840  0.3360    0.3292  0.2677   1.0000     0.5096   0.0717    0.0809
# s(Temp700)    0   0.1490  0.5654    0.6400  0.3326   0.3353     1.0000   0.1331    0.0996
# s(Slope)      0   0.0577  0.1944    0.1940  0.2844   0.0945     0.3463   1.0000    0.2540
# s(Aspect)     0   0.0194  0.1047    0.1216  0.1523   0.0382     0.1043   0.1911    1.0000


# Sal0 problematic w SSH, Sal700, Temp0, less w Chl0
# Sal700 prolematic w Sal0, Temp700
# Taking out Sal0, Sal700, SSH

# run reduced model
redMod = gam(Pres ~ s(FSLE0,bs="cs",k=5)
             # + s(Sal0,bs="cs",k=5)
             # + s(Sal700,bs="cs",k=5)
             # + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=5)
             + s(Temp700,bs="cs",k=4)
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

#             para s(FSLE0) s(Temp0) s(Temp700) s(Slope) s(Aspect)
# para          1   0.0000   0.0000     0.0000   0.0000    0.0000
# s(FSLE0)      0   1.0000   0.0448     0.2162   0.0401    0.0267
# s(Temp0)      0   0.0840   1.0000     0.5096   0.0717    0.0809
# s(Temp700)    0   0.1490   0.3353     1.0000   0.1331    0.0996
# s(Slope)      0   0.0577   0.0945     0.3463   1.0000    0.2540
# s(Aspect)     0   0.0194   0.0382     0.1043   0.1911    1.0000

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
#                                          k = 5) + s(Slope, bs = "cs", k = 5) + s(Temp0, bs = "cs", 
#                                                                                  k = 5) + s(Temp700, bs = "cs", k = 4) + 1
# 
# Parametric coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.18277    0.02782  -42.52   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                 edf Ref.df  Chi.sq p-value    
#   s(Aspect)  2.992      3  3838.1  <2e-16 ***
#   s(FSLE0)   3.694      4   254.6  <2e-16 ***
#   s(Slope)   3.994      4 22391.1  <2e-16 ***
#   s(Temp0)   3.976      4   675.6  <2e-16 ***
#   s(Temp700) 2.969      3   534.4  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.574   Deviance explained = 72.3%
# -REML =  19564  Scale est. = 1         n = 10477


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
                      # + s(Sal0,bs="cs",k=5)
                      # + s(Sal700,bs="cs",k=5)
                      # + s(SSH0,bs="cs",k=5)
                      + s(Temp0,bs="cs",k=5)
                      + s(Temp700,bs="cs",k=4),
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
    
    
      # Remove non-significant terms & re-run model iteratively until only signif covars remain
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
        sitePV = summary(redMod)$s.pv
        if (!any(redSitePV>=0.05)){
          siteModList[[sites[i]]] = optSiteMod
          pValList[[sites[i]]] = sitePV
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

