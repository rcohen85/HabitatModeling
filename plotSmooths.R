plotSmooths = function(mod,data,covar,coefInd,k,periodic,resids,site,title){
  
  probs = list(c(0.5),
            c(0.333,0.666),
            c(0.275,0.5,0.725),
            c(0.2,0.4,0.6,0.8))
  quant.func<- function(x){quantile(x, probs=c(0.0275,0.975))}
  
  if (!is.na(site)){
    ind = which(!is.na(str_match(data$Site,site)))
  } else {
    ind = 1:dim(data)[1]
  }
  
  minVal = min(data[[covar]][ind])
  maxVal = max(data[[covar]][ind])
  varSeq = seq(minVal,maxVal,length.out=1000)
  varBasis = mSpline(varSeq,
                     knots=quantile(data[[covar]][ind],probs=unlist(probs[k])),
                     Boundary.knots=c(minVal,maxVal),
                     periodic=periodic)
  
  Fit = varBasis%*%coef(mod)[coefInd]
  realFit = Fit+coef(mod)[1]
  plotDF = data.frame(Var=varSeq,Fit=realFit)
  
  bootstrapParams = rmvnorm(10000, coef(mod), summary(mod)$cov.unscaled)
  bootstrapCoefs = bootstrapParams[,coefInd]
  
  bootstrapFits = (varBasis%*%t(bootstrapCoefs))+coef(mod)[1]
  cis<-apply(bootstrapFits, 1, quant.func)
  
  if (!resids){
    ggplot(plotDF, aes(Var, Fit),
    ) + geom_smooth(aes(ymin=cis[1,], ymax=cis[2,]),
                    color="#16A7CA",
                    fill="#16A7CA",
                    alpha=0.3,
                    stat ="identity"
    ) + geom_rug(data=data[ind,],
                 inherit.aes=F,
                 aes(x=.data[[covar]][ind]),
                 sides="b"
    ) + labs(x = covar,
             y = paste('s(',covar,')',sep="")
    ) + ggtitle(title
    ) + theme(axis.line = element_line(size=0.2),
              panel.background = element_blank()
    )
    
  }else if (resids){
    ggplot(plotDF, aes(Var, Fit),
    ) + geom_point(data=data[ind,],
                   aes(x=.data[[covar]][ind],y=log10(Pres)),
                   size=0.75,
                   color="#71797E"
    ) + geom_smooth(aes(ymin=cis[1,], ymax=cis[2,]),
                    color="#16A7CA",
                    fill="#16A7CA",
                    alpha=0.3,
                    stat ="identity"
    ) + geom_rug(data=data[ind,],
                 inherit.aes=F,
                 aes(x=.data[[covar]][ind]),
                 sides="b"
    ) + labs(x = covar,
             y = paste('s(',covar,')',sep="")
    ) + ggtitle(title
    ) + theme(axis.line = element_line(size=0.2),
              panel.background = element_blank()
    )
  }
  
}