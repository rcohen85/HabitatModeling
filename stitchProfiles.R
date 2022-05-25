# NOTE: If remaking profile files, delete old ones from inDir, or they will be concatenated with new ones
library(pracma)
inDir = 'J:/Chpt_3/CovarTS'
fileList = dir(inDir)
covar = "EKE"
lon = "ES"
depths=c(seq(0,100,by=10),seq(150,300,by=50),seq(400,800,by=100))

whichInd = which(!is.na(str_match(fileList,covar)) & !is.na(str_match(fileList,lon)))
HZProfile = numeric()
OCProfile = numeric()
NCProfile = numeric()
BCProfile = numeric()
WCProfile = numeric()
NFCProfile = numeric()
HATProfile = numeric()
GSProfile = numeric()
BPProfile = numeric()
BSProfile = numeric()
JAXProfile = numeric()

for (i in 1:length(depths)){
  
  thisDepth = which(!is.na(str_match(fileList[whichInd],paste('_',as.character(depths[i]),'_',sep=""))))
  if (!isempty(thisDepth)){
  load(paste(inDir,'/',fileList[whichInd[thisDepth]],sep=""))
  } else {
      masterData.Data = rep(NA,nrow=11,ncol=1191)
    }
  
  HZProfile = rbind(HZProfile,masterData.Data[1,])
  OCProfile = rbind(OCProfile,masterData.Data[2,])
  NCProfile = rbind(NCProfile,masterData.Data[3,])
  BCProfile = rbind(BCProfile,masterData.Data[4,])
  WCProfile = rbind(WCProfile,masterData.Data[5,])
  NFCProfile = rbind(NFCProfile,masterData.Data[6,])
  HATProfile = rbind(HATProfile,masterData.Data[7,])
  GSProfile = rbind(GSProfile,masterData.Data[8,])
  BPProfile = rbind(BPProfile,masterData.Data[9,])
  BSProfile = rbind(BSProfile,masterData.Data[10,])
  JAXProfile = rbind(JAXProfile,masterData.Data[11,])
}

save(HZProfile,
     OCProfile,
     NCProfile,
     BCProfile,
     WCProfile,
     NFCProfile,
     HATProfile,
     GSProfile,
     BPProfile,
     BSProfile,
     JAXProfile,
     masterData.Time,
     depths,
     file=paste(inDir,'/',covar,'_Profiles_ES.Rdata',sep=""))
