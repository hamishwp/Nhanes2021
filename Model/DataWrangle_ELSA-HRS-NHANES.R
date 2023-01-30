load("/home/patten/Documents/Coding/Oxford/BP_NhanesA/ELSA_BloodPressure.Rdata")
ELSA<-backgroundSystDiast; rm(backgroundSystDiast)
ELSA$sys1<-as.numeric(ELSA$sys1) ; ELSA$dias1<-as.numeric(ELSA$dias1)

ELSA%<>%mutate(MeanS=rowMeans(cbind(sys1,sys2,sys3),na.rm=T),
               sdS=apply(cbind(sys1,sys2,sys3),1,sd,na.rm=T),
               MeanD=rowMeans(cbind(dias1,dias2,dias3),na.rm=T),
               sdD=apply(cbind(dias1,dias2,dias3),1,sd,na.rm=T))

tmp<-ELSA%>%group_by(idauniq)%>%summarise(diffMS=MeanS-MeanS[which.min(visyear)],
                                          diffsdS=sdS-sdS[which.min(visyear)],
                                          diffMD=MeanD-MeanD[which.min(visyear)],
                                          diffsdD=sdD-sdD[which.min(visyear)],.groups="keep")
ELSA<-cbind(ELSA,tmp[,2:5]); rm(tmp)

ggplot(ELSA)+geom_point(aes(MeanS,diffMS))
