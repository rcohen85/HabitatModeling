# Cluster covariate conditions across sites to determine sensible groupings
library(useful)
library(cluster)
library(ggbiplot)

masterDF = data.frame(read.csv('MasterWeeklyDF.csv'))
# remove NAs
badRows = unique(which(is.na(masterDF),arr.ind=TRUE)[,1])
masterDF = masterDF[-badRows,]

# create training data with no site labels
trainDat = masterDF[,3:14]

# determine optimal # clusters
HartsRule = FitKMeans(trainDat, max.clusters=20, nstart=25)
PlotHartigan(HartsRule)

# also look at gap statistic
theGap = clusGap(trainDat, FUNcluster=pam, K.max=20) 
gapDF = as.data.frame(theGap$Tab)
ggplot(gapDF, aes(x=1:nrow(gapDF))
       )+geom_line(aes(y=gap), color="red"
       )+geom_point(aes(y=gap), color="red"
       )+geom_errorbar(aes(ymin=gap-SE.sim, ymax=gap+SE.sim), color="red"
       )+labs(x="Number of Clusters", y="Gap")

save(HartsRule,theGap,file='ClusterMetrics.Rdata')

# perform k-means clustering
K = kmeans(trainDat,centers=3,nstart=25)

# plot, coding by site and cluster
plot(K, data=masterDF[,2:14], class="Site")


# just cluster median covar values
medMat = matrix(nrow=length(unique(masterDF$Site)),ncol=)

