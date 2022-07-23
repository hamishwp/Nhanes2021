library(pracma)
library(plyr)
library(magrittr)
# These commands starts the parallel workers
library(rstan)
#options(mc.cores = parallel::detectCores())
options(mc.cores = 1)

args <- commandArgs(trailingOnly = TRUE)
#
# `a_filename` is a full path like 
#  "../Samples/savedRun_nhanesFRS_basehazEst_dias_eventall_FRS1998.csv"
#  The first argument to the Rscript will be used as the value for a_filename

directory<- "/data/localhost/patten/BP_NhanesA/"
#name=paste0("savedRun_nhanesA__sys_eventCVDHrt_5k_sigma_2_",args[1])
name=args[1]

a_filename <- paste0(directory,"Samples/",name,".csv") #args[1]

# Use a regex to extract the bit between "../Samples/<I-want-this-string>.csv
filename_to_save <- name #gsub("^.*/(.*)\\.csv", "\\1", a_filename)

# Keep a log of the output
sink(file=paste0(name,".txt"))

##################### Linear Predictor Calculation ########################

##################### Linear Predictor Calculation ########################
Survival_NLL_NF<-function(M_i_S,M_i_D,D_i_S,D_i_D,tau_C_S,tau_C_D,tau_H_S,tau_H_D,
                          age,delta,Time,beta,B,theta,xhat){
  
  xhat%<>%exp()%>%array(dim=c(4,2))
  
  dimmy<-dim(D_i_S)
  
  # beta[ 1:1000 ]
  # x[ 1:1000  , 1:18018 ]
  A1<-beta[,1]*(M_i_S - xhat[1,1])
  A2<-beta[,2]*(abs(D_i_S) - xhat[2,1])
  A3<-beta[,3]*(tau_C_S - xhat[3,1])
  A4<-beta[,4]*(tau_H_S - xhat[4,1])
  A5<-beta[,5]*(M_i_D - xhat[1,2])
  A6<-beta[,6]*(abs(D_i_D) - xhat[2,2])
  A7<-beta[,7]*(tau_C_D - xhat[3,2])
  A8<-beta[,8]*(tau_H_D - xhat[4,2])  
  
  linpred<- A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8
  
  D<-mean(abs(rowSums(exp(linpred))-dimmy[2])/dimmy[2])
  
  D1<-mean(abs(rowSums(exp(A1))-dimmy[2])/dimmy[2])
  D2<-mean(abs(rowSums(exp(A2))-dimmy[2])/dimmy[2])
  D3<-mean(abs(rowSums(exp(A3))-dimmy[2])/dimmy[2])
  D4<-mean(abs(rowSums(exp(A4))-dimmy[2])/dimmy[2])
  D5<-mean(abs(rowSums(exp(A5))-dimmy[2])/dimmy[2])
  D6<-mean(abs(rowSums(exp(A6))-dimmy[2])/dimmy[2])
  D7<-mean(abs(rowSums(exp(A7))-dimmy[2])/dimmy[2])
  D8<-mean(abs(rowSums(exp(A8))-dimmy[2])/dimmy[2])
  
  rm(A1,A2,A3,A4,A5,A6,A7)
  
  funcy<-8*D+D1+D2+D3+D4+D5+D6+D7+D8
  
  full<-abs(mean(exp(linpred))-1)
  print(full)
  
  n<-20
  breaks<-seq(from=min(age,na.rm = T)-1, to=max(age,na.rm = T)+1,length.out=n)
  a_j<-H_j<-d_j<-vector(length = (n-2))
  
  for (i in 2:(n-1)){
    
    # Age histogram a_j
    a_j[i-1]<-breaks[i]
    # indices of people that lived through age a_j
    ind<-age<=breaks[i]

    # cumulative hazard of all people that lived through age a_j
    H_j[i-1]<-sum(colMeans((B[,ind]/theta[,ind])*exp(age[ind]*theta[,ind] + linpred[,ind])*(exp(theta[,ind]*Time[ind]) - 1),na.rm = T),na.rm = T)
    # cumulative deaths of all people that lived up to age a_j
    d_j[i-1]<-sum(delta[ind])
  }
  
  diffy<-mean(sqrt((H_j-d_j)*(H_j-d_j)))
  print(diffy)
  print("")
  return(diffy*funcy*full)
  
}

Survival_NLL_F<-function(FRS,FRSc,D_i_S,D_i_D,tau_C_S,tau_C_D,tau_H_S,tau_H_D,
                         age,delta,Time,beta,B,theta,xhat){

  xhat<-c(0,c(xhat)[1:3],0,c(xhat)[4:6])
  xhat%<>%exp()%>%array(dim=c(4,2))
  
  dimmy<-dim(D_i_S)
  
  # beta[ 1:1000 ]
  # x[ 1:1000  , 1:18018 ]
  A1<-beta[,1]%o%(FRS-FRSc)
  A2<-beta[,2]*(abs(D_i_S) - xhat[2,1])
  A3<-beta[,3]*(tau_C_S - xhat[3,1])
  A4<-beta[,4]*(tau_H_S - xhat[4,1])
  A5<-beta[,5]*(abs(D_i_D) - xhat[2,2])
  A6<-beta[,6]*(tau_C_D - xhat[3,2])
  A7<-beta[,7]*(tau_H_D - xhat[4,2])
  
  linpred<-A1 + A2 + A3 + A4 + A5 + A6 + A7
  
  D<-mean(abs(rowSums(exp(linpred))-dimmy[2])/dimmy[2])
  
  D2<-mean(abs(rowSums(exp(A2))-dimmy[2])/dimmy[2])
  D3<-mean(abs(rowSums(exp(A3))-dimmy[2])/dimmy[2])
  D4<-mean(abs(rowSums(exp(A4))-dimmy[2])/dimmy[2])
  D5<-mean(abs(rowSums(exp(A5))-dimmy[2])/dimmy[2])
  D6<-mean(abs(rowSums(exp(A6))-dimmy[2])/dimmy[2])
  D7<-mean(abs(rowSums(exp(A7))-dimmy[2])/dimmy[2])
  
  rm(A1,A2,A3,A4,A5,A6,A7)
  
  funcy<-6*D+D2+D3+D4+D5+D6+D7
  
  full<-abs(mean(exp(linpred))-1)
  print(full)
  
  n<-20
  breaks<-seq(from=min(age,na.rm = T)-1, to=max(age,na.rm = T)+1,length.out=n)
  a_j<-H_j<-d_j<-vector(length = (n-2))
  
  for (i in 2:(n-1)){
    
    # Age histogram a_j
    a_j[i-1]<-breaks[i]
    # indices of people that lived through age a_j
    ind<-age<=breaks[i]
    
    # cumulative hazard of all people that lived through age a_j
    H_j[i-1]<-sum(colMeans((B[,ind]/theta[,ind])*exp(age[ind]*theta[,ind] + linpred[,ind])*(exp(theta[,ind]*Time[ind]) - 1),na.rm = T),na.rm = T)
    
    # cumulative deaths of all people that lived up to age a_j
    d_j[i-1]<-sum(delta[ind])
  }
  
  diffy<-mean(sqrt((H_j-d_j)*(H_j-d_j)))
  print(diffy)
  print("")
  return(diffy*funcy*full)
  
}
####################### Bayes Factor Calculation ##########################
bayesFac <- function(mcmcobj){
  
  if(is.list(mcmcobj)) return(ldply(mcmcobj, bayesFac))
  
  if(is.array(mcmcobj) &&
     !is.na(ncol(mcmcobj)) &&
     !is.null(ncol(mcmcobj)) &&
     (ncol(mcmcobj) > 1)) return( adply(mcmcobj, 2, bayesFac) )
  
  mean_estimate <- mean(mcmcobj)
  larger <- (mean_estimate >= 0)
  p_larger_than_0 <- mean(mcmcobj>=0)
  p_smaller_than_0 <- 1 - p_larger_than_0
  
  if(larger){
    h1_string <- "par > 0"
    bayesfac <- p_larger_than_0/p_smaller_than_0
  } else {
    h1_string <- "par < 0"
    bayesfac <- p_smaller_than_0/p_larger_than_0
  }
  
  return(data.frame("H_1"=h1_string, "BF"=bayesfac))
  
}


############### Extract Results & Calc Centering Values ###################
make_results <- function(a_filename, listnames){
  
  load(listnames)
  
  if(grepl(a_filename,pattern = "FRS",fixed = T)) {
    list_nhanes<-list_nhanesFRS
  } else list_nhanes<-list_nhanesA
  
  Time<-list_nhanes$T
  age<-list_nhanes$age
  gender<-list_nhanes$male+1
  ethn<-list_nhanes$black+2*list_nhanes$white+3*list_nhanes$other
  
  FR_score<-grepl(a_filename,pattern = "ATP",fixed = T)|grepl(a_filename,pattern = "1998",fixed = T)
  
  if (FR_score){
    
    if(grepl(a_filename,pattern = "ATP",fixed = T)) {
      FRS<-list_nhanes$FRS.ATP
    } else FRS<-list_nhanes$FRS.1998
    FRSc<-log(mean(exp(FRS)))
    
    if(grepl(a_filename,pattern = "CVDHrt",fixed = T)) {
      delta<-list_nhanes$eventCVDHrt
    } else delta<-list_nhanes$eventall
    
  } else {
    
    if(grepl(a_filename,pattern = "CVDHrt",fixed = T)) {
      delta<-list_nhanes$eventCVDHrt
    } else delta<-list_nhanes$eventall
    
  }
  
  cat("\n\nTrying to read", a_filename, "...")
  fit <- read_stan_csv(a_filename)
  cat("success.\n\n")
  
  beta<-as.matrix(fit,pars=c("beta"))
  gompertz<-as.matrix(fit,pars=c("gompertz"))
  
  MMM<-as.matrix(fit,pars=c("M_i"))
  Delt<-as.matrix(fit,pars=c("Delta_i"))
  tauis<-as.matrix(fit,pars=c("tauis"))
  fit.summary <- summary(fit)
  fit.summary <- fit.summary$summary
  BFs <- bayesFac(extract(fit, pars=c("beta")))
  times <- get_elapsed_time(fit)
  rm(fit)
  
  M_i_S <-MMM[,grep("M_i.*,1]",colnames(MMM))] # Systolic clinic-home mean BP
  M_i_D <-MMM[,grep("M_i.*,2]",colnames(MMM))] # Diastolic clinic-home mean BP
  rm(MMM)
           
  D_i_S <-Delt[,grep("Delta_i.*,1]",colnames(Delt))] # Systolic clinic-home difference BP
  D_i_D <-Delt[,grep("Delta_i.*,2]",colnames(Delt))] # Diastolic clinic-home difference BP
  rm(Delt)
           
  tau_C_S <-tauis[,grep("tau.*,1]",colnames(tauis))] # Systolic clinic precision
  tau_H_S <-tauis[,grep("tau.*,2]",colnames(tauis))] # Systolic home precision
  tau_C_D <-tauis[,grep("tau.*,3]",colnames(tauis))] # Diastolic clinic precision
  tau_H_D <-tauis[,grep("tau.*,4]",colnames(tauis))] # Diastolic home precision
  rm(tauis)
  
  gompz<-array(gompertz, c(dim(gompertz)[1],2,3,2))
  
  B<-theta<-array(0,dim=dim(D_i_S))
  for (j in 1:length(gender)){
    B[,j]<-gompz[,gender[j],ethn[j],1]
    theta[,j]<-gompz[,gender[j],ethn[j],2]
  }
  rm(gompz,gender,ethn)
    
  if(!FR_score){
    
    # xhat<-array(c(128,5.3,0.075,0.115,  75,3.85,0.15,0.19),dim = c(4,2))
    xhat <-array(c(mean(M_i_S),mean(D_i_S),mean(tau_C_S),mean(tau_H_S),
                  mean(M_i_D),mean(D_i_D),mean(tau_C_D),mean(tau_H_D)),dim = c(4,2))
    xhat%<>%log()
    func <- function(xh) Survival_NLL_NF(M_i_S=M_i_S,M_i_D=M_i_D,D_i_S=D_i_S,D_i_D=D_i_D,
                                           tau_C_S=tau_C_S,tau_C_D=tau_C_D,tau_H_S=tau_H_S,tau_H_D=tau_H_D,
                                           age=age,delta=delta,Time=Time,beta=beta,B=B,theta=theta,xhat=xh)

    output<-optim(xhat,func)
    xhat<-exp(output$par)
    xhat[2,]%<>%abs
    
  }
  else{
    
    # xhat<-array(c(0,5.3,0.075,0.115,  0,3.85,0.15,0.19),dim = c(4,2))
    xhat<-array(c(FRSc,mean(D_i_S),mean(tau_C_S),mean(tau_H_S),
                  1e-6,mean(D_i_D),mean(tau_C_D),mean(tau_H_D)),dim = c(4,2))
    xhat%<>%log()
      
    func <- function(xh) Survival_NLL_F(FRS=FRS,FRSc=FRSc,D_i_S=D_i_S,D_i_D=D_i_D,
                                          tau_C_S=tau_C_S,tau_C_D=tau_C_D,tau_H_S=tau_H_S,tau_H_D=tau_H_D,
                                          age=age,delta=delta,Time=Time,beta=beta,B=B,theta=theta,xhat=xh)
    
    output<-optim(xhat[-1,],func)
    xhat[-1,]<-exp(output$par)
    xhat[2,]%<>%abs
    
  } 
  
  print("xhat values:")
  print(xhat[1,])
  print(xhat[2,])
  print(xhat[3,])
  print(xhat[4,])
  print(" ")
  
  print("NLL:")
  print(output$value)
  
  resultslist <- list()
  resultslist$fit.summary=fit.summary
  resultslist$set=a_filename
  resultslist$beta=beta
  resultslist$gompertz=gompertz
  resultslist$BF<-BFs
    
  resultslist$M_i_S<-M_i_S
  resultslist$M_i_D<-M_i_D
  resultslist$D_i_S<-D_i_S
  resultslist$D_i_D<-D_i_D
  resultslist$tau_H_S<-tau_H_S
  resultslist$tau_C_S<-tau_C_S
  resultslist$tau_H_D<-tau_H_D
  resultslist$tau_C_D<-tau_C_D
    
  resultslist$xhat1=xhat[1,]
  resultslist$xhat2=xhat[2,]
  resultslist$xhat3=xhat[3,]
  resultslist$xhat4=xhat[4,]
  
  resultslist$NLL<-output$value
  
  return(resultslist)
  
}

resultslist <- make_results(a_filename,listnames=paste0(directory,"Data_cleaned/nhanes_cleaned_lists.RData"))

save(resultslist, file=paste0(directory,"Samples/summariesdata/",filename_to_save,".Rdata"))

#Close the sink
sink()

