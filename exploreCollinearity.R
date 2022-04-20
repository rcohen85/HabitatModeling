# Procedure to compare covariates and assess collinearity prior to constructing habitat models. 
# VIF calculated based on methods and code from the usdm package

# library(usdm)
library(viridis)
library(stringr)
library(ggplot2)
saveName = 'SurfPlus600'

data = as.data.frame(read.csv('J:/Chpt_3/ModelData/Gervais_masterDF.csv'))
vars = data.frame(colnames(data))
tempDF = data[,c(4:21,64:70,78:84,127:133,183:189,239:245)]
# tempInd = which(!is.na(str_match(colnames(data),"Temp")))
# tempDF = data[,tempInd]
VIFmat = matrix(ncol=ncol(tempDF),nrow=ncol(tempDF))
colnames(VIFmat) = colnames(tempDF)
rownames(VIFmat) = colnames(tempDF)


for (i in 1:ncol(tempDF)){ # for each covar (column)
  resp = tempDF[,i]

  for (k in 1:ncol(tempDF)){ # calculate VIF relative to each other covar
    
  # create linear models
    lm1 = lm(resp~tempDF[,k])
    lm2 = lm(tempDF[,k]~resp)
    
  # calculate VIF from R2
    VIF = c(1/(1-summary(lm1)$r.squared),1/(1-summary(lm2)$r.squared))
  
  # save larger VIF value
    mn = which.max(VIF)
    VIFmat[i,k] = VIF[mn]
    
  }
}


# plot VIF values as heatmap

plotDF = data.frame(z=stack(data.frame(VIFmat))[,1],
                    x=rep(colnames(tempDF),times=ncol(tempDF)),
                    y=rep(colnames(tempDF),each=ncol(tempDF)))
plotDF$z[plotDF$z>=20] = 20

ggplot(plotDF,aes(x=x,y=y)
  ) + geom_tile(aes(fill=z)
  ) + scale_fill_viridis(
  ) + theme(axis.text.x = element_text(angle = 60,hjust=1)
  ) + labs(x=NULL,y=NULL
  ) + scale_x_discrete(limits=colnames(tempDF)
  ) + scale_y_discrete(limits=colnames(tempDF))
ggsave(filename=paste('J:/Chpt_3/ModelData/',saveName,'.png',sep=""),scale=2)



