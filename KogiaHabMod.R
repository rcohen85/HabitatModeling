library(tidyverse)
# library(mgcv.helper)
library(splines2)
library(mgcv)
library(MuMIn)
library(gratia)

## GAM approach ---------------------
# Regional model
spec = 'Kogia'
outDir = "E:/ModelingCovarData/ModelOutput"

# if it doesn't already exist, create directory to save models and figures
if (!dir.exists(paste(outDir,'/',spec,sep=""))){
  dir.create(paste(outDir,'/',spec,sep=""))
}

data = data.frame(read.csv('E:/ModelingCovarData/Master_DFs/Kogia_masterDF.csv'))
# Round presence to get Poisson dist
data$Pres = round(data$Pres)

# Check which covars are correlated w presence to determine starting covar list

# Test for how a term should be included in the model
smoothVarList = c ("Chl0",
                   "EKE0",
                   "FSLE0",
                   "Sal0",
                   "Sal700",
                   "SSH0",
                   "Temp0",
                   "Temp700",
                   "VelAsp0",
                   "VelMag0",
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
# Chl0    "14490.731614196"  "14320.7573161558" "14319.7075979385" "14063.5069285766" "fiveKnots"
# EKE0    "16823.2975907724" "16692.8474205127" "16541.7308311942" "16529.804451898"  "fiveKnots"
# FSLE0   "17041.4033417367" "16870.595490791"  "16863.6867957874" "16862.6645863488" "fiveKnots"
# Sal0    "14181.7189711522" "13847.9819294297" "13805.3806250525" "13790.1401037295" "fiveKnots"
# Sal700  "14046.3503845738" "13842.8772721015" "13826.2353758022" "13806.5799752475" "fiveKnots"
# SSH0    "13858.8509733107" "13617.9804640047" "13457.4544634775" "13444.0205139895" "fiveKnots"
# Temp0   "15500.5624432781" "15314.7646282148" "15201.5488772175" "15218.4316693565" "fourKnots"
# Temp700 "14590.7159920157" "14366.7683365063" "14078.3752650435" "13991.660634083"  "fiveKnots"
# VelAsp0 "16837.7250307821" "16553.9988199736" "16553.9988199736" "16537.8297596387" "fiveKnots"
# VelMag0 "17102.1959279262" "16856.1743787049" "16796.3808693215" "16780.9064501885" "fiveKnots"
# Slope   "14925.7394366375" "14869.2952437645" "14511.4868173416" "14273.4287770369" "fiveKnots"
# Aspect  "15636.9138553194" "16695.1980455299" "16695.1980455299" "16394.0574411581" "linMod"   

# run full model
fullMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
              + s(EKE0,bs="cs",k=5)
              + s(FSLE0,bs="cs",k=5)
              + s(Sal0,bs="cs",k=4)
              + s(Sal700,bs="cs",k=5)
              + s(SSH0,bs="cs",k=5)
              + s(Temp0,bs="cs",k=4)
              + s(Temp700,bs="cs",k=5)
              + s(VelAsp0,bs="cc",k=5)
              + s(VelMag0,bs="cs",k=5)
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

#             para s(Chl0) s(EKE0) s(FSLE0) s(Sal0) s(Sal700) s(SSH0) s(Temp0) s(Temp700) s(VelAsp0) s(VelMag0) s(Slope) s(Aspect)
# para          1  0.0000  0.0000   0.0000  0.0000    0.0000  0.0000   0.0000     0.0000     0.0000     0.0000   0.0000    0.0000
# s(Chl0)       0  1.0000  0.0776   0.0322  0.3345    0.2248  0.2111   0.4676     0.2390     0.0699     0.0366   0.1385    0.0477
# s(EKE0)       0  0.0315  1.0000   0.1023  0.2900    0.1849  0.1574   0.1456     0.1886     0.9430     0.2263   0.0278    0.0395
# s(FSLE0)      0  0.0206  0.0812   1.0000  0.1710    0.1125  0.1191   0.0817     0.1255     0.0712     0.1075   0.0401    0.0267
# s(Sal0)       0  0.1286  0.1797   0.1268  1.0000    0.6319  0.2476   0.3911     0.5615     0.1608     0.1795   0.1414    0.1053
# s(Sal700)     0  0.1230  0.2105   0.1531  0.8001    1.0000  0.3452   0.3862     0.6636     0.1918     0.2859   0.1839    0.1746
# s(SSH0)       0  0.1551  0.2027   0.1434  0.5754    0.3864  1.0000   0.4933     0.4568     0.1849     0.3112   0.2119    0.1298
# s(Temp0)      0  0.2273  0.1145   0.0835  0.3614    0.3221  0.2631   1.0000     0.3928     0.1044     0.1349   0.0633    0.0635
# s(Temp700)    0  0.1252  0.1993   0.1492  0.7265    0.6841  0.3361   0.4341     1.0000     0.1801     0.1971   0.1361    0.1256
# s(VelAsp0)    0  0.0318  0.8862   0.0981  0.2848    0.1819  0.1511   0.1447     0.1861     1.0000     0.2129   0.0236    0.0366
# s(VelMag0)    0  0.0171  0.1400   0.1150  0.2339    0.1524  0.2408   0.1108     0.1385     0.1297     1.0000   0.0208    0.0959
# s(Slope)      0  0.0919  0.0649   0.0577  0.2556    0.1940  0.2844   0.1489     0.1849     0.0607     0.0767   1.0000    0.2540
# s(Aspect)     0  0.0361  0.0348   0.0194  0.1489    0.1216  0.1523   0.0573     0.0618     0.0419     0.1689   0.1911    1.0000

# SSH problematic with Chl0
# EKE0 problematic with VelAsp0
# Sal0 problematic w Sal700, Temp700, less w SSH, Temp0
# Sal700 problematic w Sal0, Temp700
# Taking out Sal0, SSH0, Temp700, VelAsp0

# run reduced model
redMod = gam(Pres ~ s(Chl0,bs="cs",k=5)
             + s(EKE0,bs="cs",k=5)
             + s(FSLE0,bs="cs",k=5)
             # + s(Sal0,bs="cs",k=4)
             + s(Sal700,bs="cs",k=5)
             # + s(SSH0,bs="cs",k=5)
             + s(Temp0,bs="cs",k=4)
             # + s(Temp700,bs="cs",k=5)
             # + s(VelAsp0,bs="cc",k=5)
             + s(VelMag0,bs="cs",k=5)
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

# sal700 still partially explained by temp0, temp0 by chl0 but proceeding anyway

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
#                                          k = 5) + s(EKE0, bs = "cs", k = 5) + s(FSLE0, bs = "cs", 
#                                                                                 k = 5) + s(Slope, bs = "cs", k = 5) + s(Temp0, bs = "cs", 
#                                                                                                                         k = 4) + s(VelMag0, bs = "cs", k = 5) + 1
# 
# Parametric coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.25939    0.04201  -53.79   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                 edf Ref.df  Chi.sq  p-value    
#   s(Aspect)  2.9648      3 363.665  < 2e-16 ***
#   s(Chl0)    1.9624      4  19.077 1.98e-05 ***
#   s(EKE0)    0.8338      4   2.922   0.0482 *  
#   s(FSLE0)   2.5734      4  20.678 3.66e-05 ***
#   s(Slope)   3.9798      4 655.730  < 2e-16 ***
#   s(Temp0)   2.7100      3 117.742  < 2e-16 ***
#   s(VelMag0) 2.5951      4  27.768 2.82e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.0864   Deviance explained = 32.9%
# -REML = 4577.5  Scale est. = 1         n = 10477


# plot
png(filename=paste(outDir,'/',spec,'/','Kogia_allSites.png',sep=""),width=600,height=600)
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
                      + s(EKE0,bs="cs",k=5)
                      + s(FSLE0,bs="cs",k=5)
                      + s(Sal700,bs="cs",k=5)
                      + s(Temp0,bs="cs",k=4)
                      + s(VelMag0,bs="cs",k=5),
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
    png(filename=paste(outDir,'/',spec,'/','Kogia_',sites[i],'.png',sep=""),width=600,height=600)
    plot.gam(siteModList[[sites[i]]],all.terms=TRUE,rug=TRUE,scale=0,pages=1)
    while (dev.cur()>1) {dev.off()}
    
  }
}

save(siteModList,pValList,siteModCompList,file=paste(outDir,'/',spec,'/','SiteSpecificModels.Rdata',sep=""))

