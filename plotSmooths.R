plotSmooths = function(mod,covar,coefInd,k,periodic,title){
  
  probs = list(c(0.5),
            c(0.333,0.666),
            c(0.275,0.5,0.725),
            c(0.2,0.4,0.6,0.8))
  quant.func<- function(x){quantile(x, probs=c(0.0275,0.975))}
  
  minVal = min(data[[covar]])
  maxVal = max(data[[covar]])
  varSeq = seq(minVal,maxVal,length.out=1000)
  varBasis = mSpline(varSeq,
                     knots=quantile(data[[covar]],probs=unlist(probs[k-2])),
                     Boundary.knots=c(minVal,maxVal),
                     periodic=periodic)
  
  Fit = varBasis%*%coef(mod)[coefInd]
  realFit = Fit+coef(mod)[1]
  plotDF = data.frame(Var=varSeq,Fit=realFit)
  
  bootstrapParams = rmvnorm(10000, coef(mod), summary(mod)$cov.unscaled)
  bootstrapCoefs = bootstrapParams[,coefInd]
  
  bootstrapFits = varBasis%*%t(bootstrapCoefs)+coef(mod)[1]
  cis<-apply(bootstrapFits, 1, quant.func)
  
  ggplot(plotDF, aes(Var, Fit),
  ) + geom_smooth(aes(ymin=cis[1,], ymax=cis[2,]),
                  color="#16A7CA",
                  fill="#16A7CA",
                  alpha=0.2,
                  stat ="identity"
  ) + geom_rug(data=data,
               inherit.aes=F,
               aes(x=.data[[covar]]),
               sides="b"
  ) + labs(x = covar,
           # y = "Probability"
           y = paste('s(',covar,')',sep="")
  ) + ggtitle(title
  ) + theme(axis.line = element_line(size=0.2),
            panel.background = element_blank()
  )
}