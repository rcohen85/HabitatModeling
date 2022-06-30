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

# K-means clustering -----------------------
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


# PCA -------------------------------------
sites = unique(masterDF$Site)
masterDF$SiteFac = factor(masterDF$Site,levels=sites)
pc = prcomp(trainDat,center=TRUE,scale=TRUE)

g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = masterDF$SiteFac,
              ellipse = TRUE,
              circle = FALSE,
              ellipse.prob = 0.66,
              varname.size=4)
g = g + scale_color_manual(values=c("#f05a43","#db51a4","#de76f5","#9b85f2",
                                    "#56b4fc","#2fc7cc","#3dd17d","#5aba30",
                                    "#bbc22f","#d6a10f"))
g <- g + theme(legend.direction = 'vertical',
               legend.position = 'right')

png(filename=paste(getwd(),'/EcologicalNichePlots/OceanographyPCA.png',sep=""),
    height=600, width=600,units="px",res=100)
print(g)
while (dev.cur()>1) {dev.off()}
