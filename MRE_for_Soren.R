library(stringr)
library(splines2)
library(geepack)

# load .cvs file with species presence and covariate observations
# data are daily; presence data are counts per day, scaled by observation effort and classification error
data = data.frame(read.csv('J:/Chpt_3/ModelData/UD28_masterDF.csv'))

# Round presence to get back to integer counts per day (Poisson dist)
data$Pres = round(data$Pres)

# create grouping variable (based on previous autocorrelation exploration)
lagID = 40
numClust = length(data$Pres)/(lagID-1)
if (numClust<length(data$Pres)){
  clustID = rep(1:ceiling(numClust),each=lagID)
  clustID = clustID[1:length(data$Pres)]
} else {
  clustID = 1:length(data$Pres)
}
data$GroupID = clustID

# Which covariates are to be included in this model (based on previous work and data exploration)
smoothVarList = c("SSH0",
                  "Chl0",
                  "Sal0",
                  "Temp0",
                  "Slope")
knotList = list(c(0.333,0.666),
                c(0.5),
                c(0.333,0.666),
                c(0.5),
                c(0.333,0.666))
linVarList = list()
smoothNameList = character()

# Make smooth terms
for (i in 1:length(smoothVarList)){
  
  if (str_detect(smoothVarList[i],"Asp")){
    eval(parse(text=paste('S_',smoothVarList[i],'= mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=unlist(knotList[i])),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')),periodic=TRUE)',sep="")))
  } else {
    eval(parse(text=paste('S_',smoothVarList[i],'= mSpline(data$',smoothVarList[i],',knots=quantile(data$',smoothVarList[i],',probs=unlist(knotList[i])),Boundary.knots=c(min(data$',smoothVarList[i],'),max(data$',smoothVarList[i],')))',sep="")))
  }
  
  smoothNameList = c(smoothNameList,paste('S_',smoothVarList[i],sep=""))
}

# Construct model formula with smooth and linear terms
thisForm = formula(paste('Pres~',paste(c(smoothNameList,linVarList),collapse="+"),sep=""))

# Run full model
fullMod = geeglm(thisForm,
                 family=poisson,
                 data=data,
                 id=GroupID,
                 corstr="unstructured")

fullMod$geese$error