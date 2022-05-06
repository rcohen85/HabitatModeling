plotLinears = function(mod,covar,coefInd,site,title){
  
  quant.func<- function(x){quantile(x, probs=c(0.0275,0.975))}
  
  if (!is.null(site)){
    ind = which(!is.na(str_match(data$Site,site)))
  } else {
    ind = 1:dim(data)[1]
  }
  
  minVal = min(data[[covar]][ind])
  maxVal = max(data[[covar]][ind])
  varSeq = seq(minVal,maxVal,length.out=1000)
  
  Fit = varSeq*coef(mod)[coefInd]
  realFit = Fit+coef(mod)[1]
  plotDF = data.frame(Var=varSeq,Fit=realFit)
  
  bootstrapParams = rmvnorm(10000, coef(mod), summary(mod)$cov.unscaled)
  bootstrapCoefs = bootstrapParams[,coefInd]
  
  bootstrapFits = (varSeq%*%t(bootstrapCoefs))+coef(mod)[1]
  cis<-apply(bootstrapFits, 1, quant.func)
  
  ggplot(plotDF, aes(Var, Fit),
  ) + geom_smooth(aes(ymin=cis[1,], ymax=cis[2,]),
                  color="#16A7CA",
                  fill="#16A7CA",
                  alpha=0.2,
                  stat ="identity"
  ) + geom_rug(data=data[ind,],
               inherit.aes=F,
               aes(x=.data[[covar]][ind]),
               sides="b"
  ) + labs(x = covar,
           y = 'log(Presence)'
  ) + ggtitle(title
  ) + theme(axis.line = element_line(size=0.2),
            panel.background = element_blank()
  )
  
}