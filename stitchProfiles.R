# NOTE: If remaking profile files, delete old ones from inDir, or they will be concatenated with new ones

inDir = 'J:/Chpt_3/CovarTS'
fileList = dir(inDir)
covar = "VelocityMag"
lon = "ES"

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

for (i in 1:length(whichInd)){
  
  load(paste(inDir,'/',fileList[whichInd[i]],sep=""))
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
     file=paste(inDir,'/',covar,'_Profiles_ES.Rdata',sep=""))
