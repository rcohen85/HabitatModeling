library(gridExtra)
library(ggplot2)
library(dplyr)

## Clean FSLE data -------------------------
load('J:/Chpt_3/FSLE_TS.Rdata')
FSLE = data.frame(HZ=masterData.Fsle[1,],OC=masterData.Fsle[2,],NC=masterData.Fsle[3,],BC=masterData.Fsle[4,],
                     WC=masterData.Fsle[5,],NFC=masterData.Fsle[6,],HAT=masterData.Fsle[7,],GS=masterData.Fsle[8,],
                     BP=masterData.Fsle[9,],BS=masterData.Fsle[10,],JAX=masterData.Fsle[11,])

sites = colnames(FSLE)

# Plot histograms of data
for (i in 1:11){
  eval(parse(text=(paste(sites[i],' = ggplot(data=FSLE)+geom_histogram(aes(x=',sites[i],'))+labs(title=sites[i])',sep=""))))
}
grid.arrange(HZ,OC,NC,BC,WC,NFC,HAT,GS,BP,BS,JAX, ncol=4,nrow=3,top="FSLE")

# Remove  extreme values (<-1) (make NA, interpolate later)
for (i in 1:11){
which_outliers = which(FSLE[,i]<(-1))
FSLE[which_outliers,i] = NA}

# Interpolate missing values
for (i in 1:11){
  skippedBins = which(is.na(FSLE[,i])) # missing data
  missFSLE = apply(cbind(FSLE[skippedBins-1,i],FSLE[skippedBins+1,i]),MARGIN=1,mean)
  FSLE[skippedBins,i] = missFSLE
}

# check for missing dates
timeDiff = diff(as.numeric(masterData.Time))
any(timeDiff>1)

# Make time lagged vectors
startInd = which(masterData.Time==as.Date('2016-05-01',origin="1970-01-01"))
for (i in 1:11){
  for (k in 1:3){
    if (k==1){
      lag = 7
      lagInd = startInd-lag
    } else if (k==2){
      lag=14
      lagInd = startInd-lag
    } else if (k==3){
      lag=21
      lagInd = startInd-lag
    }
  eval(parse(text=paste('FSLE$',sites[i],'Lag',lag,' = NA',sep="")))
  eval(parse(text=paste('FSLE$',sites[i],'Lag',lag,'[startInd:1216] = FSLE$',sites[i],'[lagInd:(1216-lag)]',sep="")))
  }
}


## Clean HYCOM data --------------------------------
library(R.matlab)

HYCOM = readMat('J:/Chpt_3/HYCOM_VerticalProfiles_forR.mat')

SSH = data.frame(HZ=t(HYCOM$SSH[[1]][[1]]),OC=t(HYCOM$SSH[[2]][[1]]),NC=t(HYCOM$SSH[[3]][[1]]),
                 BC=t(HYCOM$SSH[[4]][[1]]),WC=t(HYCOM$SSH[[5]][[1]]),NFC=t(HYCOM$SSH[[6]][[1]]),
                 HAT=t(HYCOM$SSH[[7]][[1]]),GS=t(HYCOM$SSH[[8]][[1]]),BP=t(HYCOM$SSH[[9]][[1]]),
                 BS=t(HYCOM$SSH[[10]][[1]]),JAX=t(HYCOM$SSH[[11]][[1]]))

vars = c('SSH','Salinity','Temperature','WaterU','WaterV')

# Histograms
for (k in vars){
  
  if (vars[k]==SSH){
    depth = 1
  } else {
    depth = c(1,15,23,28,33) # desired depth layers
  }
  
  for (j in depth){
    for (i in 1:11){
      
      eval(parse(text=(paste(sites[i],' = ggplot(data=FSLE)+geom_histogram(aes(x=',sites[i],'))+labs(title=sites[i])',sep=""))))
      
      if (i==11){
        grid.arrange(HZ,OC,NC,BC,WC,NFC,HAT,GS,BP,BS,JAX, ncol=4,nrow=3,top=paste())
      }
    }
  }
}


# Remove outliers

# Interpolate missing values

# Create time lagged vectors


## Make modeling csv files ---------------------------

