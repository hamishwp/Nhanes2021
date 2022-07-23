source("./Model/Functions.R")

# Extract xhat values, all beta, gompertz and BP data, and save as iteration + 1
folder<-"/home/patten/Documents/Coding/Oxford/BP_NhanesA/Samples/Tau_autopred_it5/"

filez<-c("savedRun_nhanesA_eventCVDHrt_tau_autopred_",
         "savedRun_nhanesA_eventall_tau_autopred_",
         "savedRun_nhanesFRSpop_eventCVDHrt_tau_autopred_",
         "savedRun_nhanesFRSpop_eventall_tau_autopred_",
         "savedRun_nhanesFRS_eventCVDHrt_tau_FRSATP_autopred_",
         "savedRun_nhanesFRS_eventall_tau_FRSATP_autopred_",
         "savedRun_nhanesFRS_eventCVDHrt_tau_FRS1998_autopred_",
         "savedRun_nhanesFRS_eventall_tau_FRS1998_autopred_")
         
# Extract and prepare RL
sigma<-F
if(sigma) variancer<-"sigma" else variancer<-"tau"
mcs<-4
it<-5
  
for(fff in filez){
  
  print(fff)
  winner<-Inf
  
  for (kk in 1:mcs){
    
    load(paste0(folder,fff,it,"_",kk,".Rdata"))
    
    RL<-resultslist
    RL$FRS<-grepl(x = fff,pattern = "FRS")
    RL$eventall<-grepl(x = fff,pattern = "eventall")
    RL$pop<-grepl(x = fff,pattern = "FRSpop")
    RL$FRSt<-NA
    if(RL$FRS & grepl(x = fff,pattern = "ATP")) RL$FRSt<-"ATP"
    else if(RL$FRS & grepl(x = fff,pattern = "1998")) RL$FRSt<-"1998"
    RL<-AddSum(RL)
    RL$fit.summary<-NULL
  
    validy<-GetSurvival(RL,roc=F,usexhat = F)
    
    if(validy$win<winner) {
      
      vsave<-validy
      wRL<-resultslist
      winner<-validy$win
      print(paste0("Winner: ",kk," with NLL=",winner))
      
    }
    
  }
  
  resultslist<-wRL; rm(wRL)
  iwin<-vsave$iwin
  
  resultslist$xhat1<-validy$xhat[1,]
  resultslist$xhat2<-validy$xhat[2,]
  resultslist$xhat3<-validy$xhat[3,]
  resultslist$xhat4<-validy$xhat[4,]
  
  resultslist$beta<-rbind(resultslist$beta[iwin,],resultslist$beta[iwin,],resultslist$beta[iwin,])
  resultslist$gompertz<-rbind(resultslist$gompertz[iwin,],resultslist$gompertz[iwin,],resultslist$gompertz[iwin,])

  resultslist$M_i_S<-rbind(resultslist$M_i_S[iwin,],resultslist$M_i_S[iwin,],resultslist$M_i_S[iwin,])
  resultslist$D_i_S<-rbind(resultslist$D_i_S[iwin,],resultslist$D_i_S[iwin,],resultslist$D_i_S[iwin,])
  resultslist$tau_C_S<-rbind(resultslist$tau_C_S[iwin,],resultslist$tau_C_S[iwin,],resultslist$tau_C_S[iwin,])
  resultslist$tau_H_S<-rbind(resultslist$tau_H_S[iwin,],resultslist$tau_H_S[iwin,],resultslist$tau_H_S[iwin,])
  
  resultslist$M_i_D<-rbind(resultslist$M_i_D[iwin,],resultslist$M_i_D[iwin,],resultslist$M_i_D[iwin,])
  resultslist$D_i_D<-rbind(resultslist$D_i_D[iwin,],resultslist$D_i_D[iwin,],resultslist$D_i_D[iwin,])
  resultslist$tau_C_D<-rbind(resultslist$tau_C_D[iwin,],resultslist$tau_C_D[iwin,],resultslist$tau_C_D[iwin,])
  resultslist$tau_H_D<-rbind(resultslist$tau_H_D[iwin,],resultslist$tau_H_D[iwin,],resultslist$tau_H_D[iwin,])
  
  resultslist$NLL<-iwin
  
  if(RL$FRS){
    if(is.na(RL$FRSt)) run_name_string<-"nhanesFRSpop"
    else run_name_string<-"nhanesFRS"
  } else run_name_string<-"nhanesA"
  
  outfile<-paste0(folder, "exp_betas", run_name_string,"_", ifelse(RL$eventall,"eventall","eventCVDHrt"),
                  ifelse(is.na(RL$FRSt), "", paste0("_FRS",RL$FRSt)),
         ifelse(sigma,"_sigma_","_tau_"),(it+1),'.Rdata')
  
  save(resultslist,file = outfile)
  
}

