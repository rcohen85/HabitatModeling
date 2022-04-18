# Procedure to compare covariates and assess collinearity prior to constructing habitat models. 
# Based on methods and code from the usdm package

# library(usdm)
library(viridis)

data = as.data.frame(read.csv('J:/Chpt_3/ModelData/Gervais_masterDF.csv'))
vars = data.frame(colnames(data))
deepDF = data[,c(4:11,40:47,76:79,108:111,140:143)]
VIFmat = matrix(ncol=ncol(deepDF),nrow=ncol(deepDF))
colnames(VIFmat) = colnames(deepDF)
rownames(VIFmat) = colnames(deepDF)


for (i in 1:ncol(deepDF)){ # for each covar (column)
  resp = deepDF[,i]

  for (k in 1:ncol(deepDF)){ # calculate VIF relative to each other covar
    
  # create linear models
    lm1 = lm(resp~deepDF[,k])
    lm2 = lm(deepDF[,k]~resp)
    
  # calculate VIF from R2
    VIF = c(1/(1-summary(lm1)$r.squared),1/(1-summary(lm2)$r.squared))
  
  # save larger VIF value
    mn = which.max(VIF)
    VIFmat[i,k] = VIF[mn]
    
  }
}


# plot VIF values as heatmap

plotDF = data.frame(z=stack(data.frame(VIFmat))[,1],
                    x=rep(colnames(deepDF),times=ncol(deepDF)),
                    y=rep(colnames(deepDF),each=ncol(deepDF)))

ggplot(plotDF,aes(x=x,y=y)
  )+geom_tile(aes(fill=z)
  )+scale_fill_viridis(
  )+theme(axis.text.x = element_text(angle = 60,hjust=1))



# usdm functions
# .vif <- function(y) {
#   z<-rep(NA,ncol(y))
#   names(z) <- colnames(y)
#   for (i in 1:ncol(y)) {
#     z[i] <-  1/(1-summary(lm(y[,i]~.,data=y[-i]))$r.squared)
#   }
#   return(z)
# }
# .vif2 <- function(y,w) {
#   z<-rep(NA,length(w))
#   names(z) <- colnames(y)[w]
#   for (i in 1:length(w)) {
#     z[i] <-  1/(1-summary(lm(as.formula(paste(colnames(y)[w[i]],"~.",sep='')),data=y))$r.squared)
#   }
#   return(z)
# }
# 
# .maxCor <- function(k){
#   k <- abs(k)
#   n <- nrow(k)
#   for (i in 1:n) k[i:n,i] <- NA
#   w <- which.max(k)
#   c(rownames(k)[((w%/%nrow(k))+1)],colnames(k)[w%%nrow(k)])
# }
# 
# 
# vifcor = function(x, th=0.9) {
#   LOOP <- TRUE
#   x <- na.omit(x)
#   n <- new("VIF")
#   n@variables <- colnames(x)
#   exc <- c()
#   while (LOOP) {
#     xcor <- abs(cor(x))
#     mx <- .maxCor(xcor)
#     if (xcor[mx[1],mx[2]] >= th) {
#       w1 <- which(colnames(xcor) == mx[1])
#       w2 <- which(rownames(xcor) == mx[2])
#       v <- .vif2(x,c(w1,w2))
#       ex <- mx[which.max(v[mx])]
#       exc <- c(exc,ex)
#       x <- x[,-which(colnames(x) == ex)]
#     } else LOOP <- FALSE
#   }
#   if (length(exc) > 0) n@excluded <- exc
#   v <- .vif(x)
#   n@corMatrix <- cor(x)
#   n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
#   n
# }
