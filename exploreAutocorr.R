library(stringr) 
library(ggfortify)
library(splines)

inDir = 'J:/Chpt_3/ModelData'
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS')

fileList = dir(inDir,pattern=".csv",full.names=TRUE)
nspec = length(fileList)

## Calculate autocorrelation of presence timeseries directly -------------------

autocorrMat = matrix(ncol=nspec,nrow=length(sites))
rownames(autocorrMat) = sites
allSpec = character()

for (i in 1:nspec){
  
  data = data.frame(read.csv(fileList[i]))
  spec = str_remove(str_remove(fileList[i],"_masterDF.csv"),paste(inDir,'/',sep=""))
  allSpec = c(allSpec,spec)
  
  for (j in 1:length(sites)){
    
    siteInd = which(!is.na(str_match(data$Site,sites[j])))
    if (sum(which(data$Pres[siteInd]>0))>10){
      acorr = acf(data$Pres[siteInd],lag.max=180,na.action=na.exclude,plot=TRUE,main=paste(spec,'at',sites[j]))
      CI = ggfortify:::confint.acf(acorr)
      ACFidx = which(acorr[["acf"]] < CI, arr.ind=TRUE)
      autocorrMat[j,i] = ACFidx[1]
      # lagID = which(abs(acorr$acf)<0.2)
      # autocorrMat[j,i] = lagID[1]
    } 
  }
}

colnames(autocorrMat) = allSpec


## Calculate autocorrelation of residuals in simple model including all possible covars ---------------


residAutocorr = matrix(ncol=nspec,nrow=1)
colnames(residAutocorr) = allSpec
allSpec = character()

for (i in 1:nspec){
  
  data = data.frame(read.csv(fileList[i]))
  spec = str_remove(str_remove(fileList[i],"_masterDF.csv"),paste(inDir,'/',sep=""))
  
  # convert hours of presence back to # 5-min bins, then round to get Poisson dist
  data$Pres = round(data$Pres*12)
  
  
  BlockMod = glm(Pres~bs(Temp0)
                 + bs(Sal0)
                 + bs(Chl0)
                 + bs(FSLE0)
                 + bs(VelMag0)
                 + bs(VelAsp0)
                 + bs(SSH0)
                 + bs(GSLat)
                 + bs(GSDist)
                 + bs(Slope)
                 + bs(Aspect),
                 data=data,family=poisson)
  
  acorr = acf(residuals(BlockMod), lag.max = 1000, main=spec)
  CI = ggfortify:::confint.acf(acorr)
  ACFidx = which(acorr[["acf"]] < CI, arr.ind=TRUE)
  residAutocorr[1,i] = ACFidx[1]
  # lagID = which(abs(acorr$acf)<0.2)
  # residAutocorr[1,i] = lagID[1]
  
}





