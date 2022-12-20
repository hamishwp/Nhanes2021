# # load("./Samples/summariesdata/savedRun_nhanesA_eventCVDHrt_sigma_0_1.Rdata")
# load("./Samples/summariesdata/savedRun_nhanesFRSpop_eventall_sigma_0_1.Rdata")

taufolder<-"/home/patten/Documents/Coding/Oxford/BP_NhanesA/Samples/Tau_manpred_it7/"
taufilers<-sort(list.files(taufolder,pattern = "exp"))
sigfolder<-"/home/patten/Documents/Coding/Oxford/Nhanes2021/Samples/summariesdata/"
sigfilers<-sort(list.files(sigfolder,pattern = "exp"))

for (i in 1:length(taufilers)){
  
  print(paste0(sigfolder,sigfilers[i]))
  load(paste0(sigfolder,sigfilers[i]))
  
  if(length(resultslist$sigma_C_S)==0){
    
    resultslist$sigma_C_S<-1./sqrt(resultslist$tau_C_S)
    resultslist$sigma_H_S<-1./sqrt(resultslist$tau_H_S)
    resultslist$sigma_C_D<-1./sqrt(resultslist$tau_C_D)
    resultslist$sigma_H_D<-1./sqrt(resultslist$tau_H_D)
    
    resultslist$tau_C_S<-NULL
    resultslist$tau_H_S<-NULL
    resultslist$tau_C_D<-NULL
    resultslist$tau_H_D<-NULL
    
  }
  
  RL<-resultslist
  
  print(paste0(taufolder,taufilers[i]))
  load(paste0(taufolder,taufilers[i]))
  
  if(ncol(resultslist$beta)==7) {
    
    RL$beta[,1]<-rep(mean(resultslist$beta[,1]),3)
    RL$beta[,2]<-rep(mean(resultslist$beta[,2]),3)
    RL$beta[,5]<-rep(mean(resultslist$beta[,5]),3)
    RL$gompertz<-resultslist$gompertz
    RL$xhat1<-resultslist$xhat1
    RL$xhat2<-resultslist$xhat2
    RL$xhat3<-c(fxhat(median(RL$sigma_C_S),mean(RL$beta[,3])),
                fxhat(median(RL$sigma_H_S),mean(RL$beta[,4])))
    RL$xhat4<-c(fxhat(median(RL$sigma_C_D),mean(RL$beta[,6])),
                fxhat(median(RL$sigma_H_D),mean(RL$beta[,7])))
    
    
  } else {
   
    RL$beta[,1]<-rep(mean(resultslist$beta[,1]),3)
    RL$beta[,2]<-rep(mean(resultslist$beta[,2]),3)
    RL$beta[,5]<-rep(mean(resultslist$beta[,5]),3)
    RL$beta[,6]<-rep(mean(resultslist$beta[,6]),3)
    RL$gompertz<-resultslist$gompertz
    RL$xhat1<-resultslist$xhat1
    RL$xhat2<-resultslist$xhat2
    RL$xhat3<-c(fxhat(median(RL$sigma_C_S),mean(RL$beta[,3])),
                fxhat(median(RL$sigma_H_S),mean(RL$beta[,4])))
    RL$xhat4<-c(fxhat(median(RL$sigma_C_D),mean(RL$beta[,7])),
                fxhat(median(RL$sigma_H_D),mean(RL$beta[,8])))
    
  }
  
  print(t(cbind(RL$xhat1,RL$xhat2,RL$xhat3,RL$xhat4)))
  
  resultslist<-RL
  save(resultslist,file = paste0(sigfolder,"mod_",sigfilers[i]))
  
}


R<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL<-resultslist
RL$FRS<-FALSE
RL$eventall<-TRUE
RL$pop<-FALSE
RL$FRSt<-NA
RL<-AddSum(RL)
RL$fit.summary<-NULL

gompz<-ConvertGompertz(RL$gompertz)
p<-ggplot(data=gompz,aes(x=Ethnicity,y=B,fill=Gender))+geom_violin()+scale_y_log10()+#+ggtitle()
  theme(plot.title = element_text(hjust = 0.5))+ylab("B - Gompertz Parameter") + xlab("Ethnicity");p
q<-ggplot(data=gompz,aes(x=Ethnicity,y=theta,fill=Gender))+geom_violin()+#+ggtitle()
  theme(plot.title = element_text(hjust = 0.5))+ylab(TeX("$\\theta$ - Gompertz Parameter")) + xlab("Ethnicity");q


RL$sigma_C_S<-1./sqrt(RL$tau_C_S)
RL$sigma_H_S<-1./sqrt(RL$tau_H_S)
RL$sigma_C_D<-1./sqrt(RL$tau_C_D)
RL$sigma_H_D<-1./sqrt(RL$tau_H_D)

if(is.na(RL$FRSt)) btmp<-ConvertBeta(RL) else {
  if(RL$FRSt=="ATP") FRS<-list_nhanesFRS$FRS.ATP else FRS<-list_nhanesFRS$FRS.1998
  btmp<-ConvertBeta(RL,FRS=FRS)
}
k<-ggplot(data=btmp,aes(x=variable,y=value,fill=variable))+geom_violin()+
  theme(plot.title = element_text(hjust = 0.5))+ylab(TeX("$\\beta$ Linear Predictor Parameter")) +
  xlab("Variable") + ggtitle(TeX("Normalised $\\beta$ Parameter"));k

