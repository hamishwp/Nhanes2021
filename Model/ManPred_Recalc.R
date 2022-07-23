source("./Model/Functions.R")

# Extract xhat values, all beta, gompertz and BP data, and save as iteration + 1
ofolder<-"/home/patten/Documents/Coding/Oxford/BP_NhanesA/Samples/Tau_autopred_it5/"
folder<-"/home/patten/Documents/Coding/Oxford/BP_NhanesA/Samples/Tau_manpred_it7/"

filez<-c("savedRun_nhanesA_eventCVDHrt_tau_",
         "savedRun_nhanesA_eventall_tau_",
         "savedRun_nhanesFRSpop_eventCVDHrt_tau_",
         "savedRun_nhanesFRSpop_eventall_tau_",
         "savedRun_nhanesFRS_eventCVDHrt_tau_FRSATP_",
         "savedRun_nhanesFRS_eventall_tau_FRSATP_",
         "savedRun_nhanesFRS_eventCVDHrt_tau_FRS1998_",
         "savedRun_nhanesFRS_eventall_tau_FRS1998_")
         
# Extract and prepare RL
sigma<-F
if(sigma) variancer<-"sigma" else variancer<-"tau"
mcs<-4
it<-7
  
for(fff in filez){
  
  print(fff)
  winner<-Inf
  
  FRS<-grepl(x = fff,pattern = "FRS")
  eventall<-grepl(x = fff,pattern = "eventall")
  pop<-grepl(x = fff,pattern = "FRSpop")
  FRSt<-NA
  if(FRS & grepl(x = fff,pattern = "ATP")) FRSt<-"ATP"
  else if(FRS & grepl(x = fff,pattern = "1998")) FRSt<-"1998"
  
  if(FRS){
    if(is.na(FRSt)) run_name_string<-"nhanesFRSpop"
    else run_name_string<-"nhanesFRS"
  } else run_name_string<-"nhanesA"
  
  oldout<-paste0(ofolder, "exp_betas", run_name_string,"_", ifelse(eventall,"eventall","eventCVDHrt"),
                 ifelse(is.na(FRSt), "", paste0("_FRS",FRSt)),
                 ifelse(sigma,"_sigma_","_tau_"),(it-1),'.Rdata')
  load(oldout)
  
  xhat1<-resultslist$xhat1
  xhat2<-resultslist$xhat2
  xhat3<-resultslist$xhat3
  xhat4<-resultslist$xhat4
  
  for (kk in 1:mcs){
    
    load(paste0(folder,fff,it,"_",kk,".Rdata"))
    
    resultslist$xhat1<-xhat1
    resultslist$xhat2<-xhat2
    resultslist$xhat3<-xhat3
    resultslist$xhat4<-xhat4
    
    RL<-resultslist
    RL$FRS<-FRS
    RL$eventall<-eventall
    RL$pop<-pop
    RL$FRSt<-FRSt
    
    RL<-AddSum(RL)
    RL$fit.summary<-NULL
  
    validy<-GetSurvival(RL,roc=F,usexhat = T)
    
    if(validy$win<winner) {
      
      vsave<-validy
      wRL<-resultslist
      winner<-validy$win
      print(paste0("Winner: ",kk," with NLL=",winner))
      
    }
    
  }
  
  plot(validy$survy$H_j,validy$survy$d_j); abline(0,1)

  resultslist<-wRL; rm(wRL)
  iwin<-vsave$iwin

  resultslist$xhat1<-xhat1
  resultslist$xhat2<-xhat2
  resultslist$xhat3<-xhat3
  resultslist$xhat4<-xhat4
  
  resultslist$NLL<-iwin
  
  outfile<-paste0(folder, "exp_betas", run_name_string,"_", ifelse(RL$eventall,"eventall","eventCVDHrt"),
                  ifelse(is.na(RL$FRSt), "", paste0("_FRS",RL$FRSt)),
         ifelse(sigma,"_sigma_","_tau_"),(it+1),'.Rdata')
  
  save(resultslist,file = outfile)
  
}

