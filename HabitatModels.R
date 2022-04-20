library(geepack)

data = data.frame(read.csv('J:/Chpt_3/ModelData/Blainville_masterDF.csv'))
data = data[order(data$Date),]
rownames(data) = seq(1,nrow(data))
ACF = acf(data$Pres,lag.max=5000)

# create grouping variable
lagID = 1300 
numClust = length(data$Pres)/(lagID-1)
if (numClust<length(data$Pres)){
  clustID = rep(1:ceiling(numClust),each=lagID)
  clustID = clustID[1:numel(data$Pres)]
} else {
  clustID = 1:length(data$Pres)
}

data$GroupID = clustID

habMod = geeglm(Pres~Temp700
                + Chl0Lag28
                + Sal700
                + FSLE0
                + SSH0
                + Slope
                + GSLat
                + VelMag700
                + Aspect:VelAsp0,
                data=data,
                family=poisson,
                id=GroupID,
                corstr="ar1")

PV = getPvalues(habMod)
