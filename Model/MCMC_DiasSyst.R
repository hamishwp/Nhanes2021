#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@            Blood Pressure Hyperparameter MCMC Calculator          @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@                            Hamish Patten                          @@@@@#
#@@@@@          Postdoctoral Researcher, University of Oxford, 2019      @@@@@#
#@@@@@                      hamish.patten@stats.ox.ac.uk                 @@@@@#
#@@@@@ Code adapted from D.B. Bester thesis, University of Oxford (2014) @@@@@#
#@@@@@            For details on code, contact Hamish Patten             @@@@@#
#@@@@@        or try David Steinsaltz at steinsal@stats.ox.ac.uk         @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

###############################################################################
######                    CODE DESCRIPTION & LAYOUT                      ######
######  - Reads in the blood pressure data of all participants           ######
######  - Reads in the MLE parameter estimates (see mortality_fit.R)     ######
######  - Launches the rstan Bayesian MCMC parameter simulation          ######
######  - Iterates on the centering values 'xhat' in log-likelihood      ######
###############################################################################

#########################        TO RUN CODE:        ##########################
# R CMD BATCH --vanilla "--args 1" MCMC_DiasSyst.R MCMC_DiasSyst.Rout &

######################### Description of Functions: ###########################
# - RunNhanesModel : given a dataset/subset, launch rstan
# - bayesFac : calculate the Bayes factor of the MCMC model
# - Fxhat : function that defines the centering parameter xhat
# - dFxhat : derivative of Fxhat
# - convert_param1 : restructure the parameter MLEs from mortality_fit.R

library(pracma)
library(plyr)
library(magrittr)
library(rstan)
rstan_options(auto_write = TRUE)

# directory<- "/data/localhost/patten/BP_NhanesA/"
directory<- "~/Documents/Coding/Oxford/BP_NhanesA/"

# Read the run number
# args <- commandArgs(trailingOnly = TRUE)
args<-c(5)

runnum <- as.integer(args[1])
print(paste0("Calculating for run number: ",runnum))

iii<-0
mcs<-4
# Set the starting point for iterative simulations of the xhat parameter
if(length(args)>1) iii<- as.integer(args[2])
print(paste0("xhat iteration number: ",iii))
# Set the number of parallel computational tasks
if(length(args)>2) mcs<- as.integer(args[3])
print(paste0("Number of parallel cpu tasks: ",mcs))

##################### Linear Predictor Calculation ########################
Survival_NLL_NF<-function(M_i_S,M_i_D,D_i_S,D_i_D,tau_C_S,tau_C_D,tau_H_S,tau_H_D,
                   age,delta,T,beta,B,theta,xhat){
  
  xhat%<>%array(dim=c(4,2))
  
  dimmy<-dim(D_i_S)
  
  # beta[ 1:1000 ]
  # x[ 1:1000  , 1:18018 ]
  A1<-beta[,1]*(M_i_S - xhat[1,1])
  A2<-beta[,2]*(abs(D_i_S) - abs(xhat[2,1]))
  A3<-beta[,3]*(tau_C_S - xhat[3,1])
  A4<-beta[,4]*(tau_H_S - xhat[4,1])
  A5<-beta[,5]*(M_i_D - xhat[1,2])
  A6<-beta[,6]*(abs(D_i_D) - abs(xhat[2,2]))
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
  breaks<-histss(age,n)$breaks
  a_j<-H_j<-d_j<-vector(length = (n-2))
  
  for (i in 2:(n-1)){
    
    # Age histogram a_j
    a_j[i-1]<-breaks[i]
    # indices of people that lived through age a_j
    ind<-age<=breaks[i]
    
    # cumulative hazard of all people that lived through age a_j
    H_j[i-1]<-sum(colMeans((B[,ind]/theta[,ind])*exp(age[ind]*theta[,ind] + linpred[,ind])*(exp(theta[,ind]*T[ind]) - 1),na.rm = T),na.rm = T)
    
    # cumulative deaths of all people that lived up to age a_j
    d_j[i-1]<-sum(delta[ind])
  }
  
  diffy<-mean(sqrt((H_j-d_j)*(H_j-d_j)))
  print(diffy)
  print("")
  return(diffy*funcy*full)
    
}

Survival_NLL_F<-function(FRS,FRSc,D_i_S,D_i_D,tau_C_S,tau_C_D,tau_H_S,tau_H_D,
                  age,delta,T,beta,B,theta,xhat){
  
  xhat<-c(0,c(xhat)[1:3],0,c(xhat)[4:6])
  xhat%<>%array(dim=c(4,2))
  
  dimmy<-dim(D_i_S)
  
  # beta[ 1:1000 ]
  # x[ 1:1000  , 1:18018 ]
  A1<-beta[,1]%o%(FRS-FRSc)
  A2<-beta[,2]*(abs(D_i_S) - abs(xhat[2,1]))
  A3<-beta[,3]*(tau_C_S - xhat[3,1])
  A4<-beta[,4]*(tau_H_S - xhat[4,1])
  A5<-beta[,5]*(abs(D_i_D) - abs(xhat[2,2]))
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
  breaks<-histss(age,n)$breaks
  a_j<-H_j<-d_j<-vector(length = (n-2))
  
  for (i in 2:(n-1)){
    
    # Age histogram a_j
    a_j[i-1]<-breaks[i]
    # indices of people that lived through age a_j
    ind<-age<=breaks[i]
    
    # cumulative hazard of all people that lived through age a_j
    H_j[i-1]<-sum(colMeans((B[,ind]/theta[,ind])*exp(age[ind]*theta[,ind] + linpred[,ind])*(exp(theta[,ind]*T[ind]) - 1),na.rm = T),na.rm = T)
    
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

set.seed(1)

############################################################################################
############################### MCMC Calculation Function ##################################
############################################################################################
RunNhanesModel <- function(my_file_stan, 
                           cleaned_nhanes_list,  # The cleaned list
                           Sestimates,
                           Destimates,
                           run_name_string,       # The string used for filename to save
                           mortality, 
                           FR_score=NULL,
                           iii=1,
                           mcs=4,
                           savelocation=paste0(directory,"Samples/")){
  ### my_file_stan - rstan code (e.g. mystanmodel_DS.stan)
  ### cleaned_nhanes_list - cleaned blood pressure measurements list (e.g. list_nhanesA)
  ### Sestimates - Systolic blood pressure hyperparameter estimates (from mortality_fit.Rmd)
  ### Destimates - Diastolic blood pressure hyperparameter estimates (from mortality_fit.Rmd)
  ### run_name_string - string for filename to save
  ### mortality - data on age and time after BP measurements were made, and whether the person died
  ### FR_score=NULL - using FRS or the full nhanesA dataset?
  ### savelocation - where to save the individual MCMC chain simulations
  ### mcs - how many Monte Carlo chains do you want to use? Gelman et al, B.D.A (2014) suggest at least 4 chains
  
  # Read in BP measurement data - note that length is number of measurements per individual (nj by N found below)
  xD <- cleaned_nhanes_list[["dias"]] # Diastolic BP measurements taken at the clinic
  xhomeD <- cleaned_nhanes_list[["dias_home"]] # Diastolic BP measurements taken at patients home
  xS <- cleaned_nhanes_list[["sys"]] # Systolic BP measurements taken at the clinic
  xhomeS <- cleaned_nhanes_list[["sys_home"]]  # Systolic BP measurements taken at patients home
  
  thisevent <- cleaned_nhanes_list[[mortality]] # Did the individual die?
  T <- cleaned_nhanes_list[["T"]] # Time since the measurements taken (has an upper limit of second census)
  nj <- cleaned_nhanes_list[["nj"]] # Number of BP measurements at the first census
  age <- cleaned_nhanes_list[["age"]] # Age at first census
  N <- cleaned_nhanes_list[["N"]] # Number of individuals in the data
  
  # FRS data set parameters - ATP or 1998 are different methods to calculate the FRS value
  if(!is.null(FR_score)){
    if(FR_score == "ATP"){
      FRS <- cleaned_nhanes_list[["FRS.ATP"]]
    } else if(FR_score == "1998") {
      FRS <- cleaned_nhanes_list[["FRS.1998"]]
    } else {
      stop("FRS score must be FRS.ATP or FRS.1998, NOTE : if using runum = 3 or 4, FR_score should be set as FALSE")
    }
    FRSc<-log(mean(exp(FRS)))
  }
  
  # Skin colour of individuals
  black <- cleaned_nhanes_list[["black"]]
  white <- cleaned_nhanes_list[["white"]]
  other <- cleaned_nhanes_list[["other"]]
  racematrix <- t(t(data.frame(black, white, other) ) * c(1L,2L,3L))
  # 1-black
  # 2-white
  # 3-other
  race <- apply(racematrix, 1, max)
  
  # Sex of individuals
  # 1-female
  # 2-male
  male <- cleaned_nhanes_list[["male"]]
  sex <- as.integer(male) + 1L
  
  # For longitudinal data analysis, we use the time-to-event
  Surv_T_and_delta <- as.array(as.matrix(data.frame(T=T, delta=thisevent)))
  
  # Initialise values
  m_M<-tau_M<-m_Delta<-tau_Delta<-alpha_CLINIC<-beta_CLINIC<-alpha_HOME<-beta_HOME<-c(0.,0.)
  
  # Hyperparameter MLE's (empirical Bayes)
  m_M[1]<-Sestimates[[1]]["m_M"]
  tau_M[1]<-1/Sestimates[[1]]["sigma2_M"]
  m_Delta[1]<-Sestimates[[1]]["m_Delta"]
  tau_Delta[1]<-1/Sestimates[[1]]["sigma2_Delta"]
  alpha_CLINIC[1]<-Sestimates[[2]]["alpha",1]
  beta_CLINIC[1]<-Sestimates[[2]]["beta",1]
  alpha_HOME[1]<-Sestimates[[2]]["alpha",2]
  beta_HOME[1]<-Sestimates[[2]]["beta",2]
  
  m_M[2]<-Destimates[[1]]["m_M"]
  tau_M[2]<-1/Destimates[[1]]["sigma2_M"]
  m_Delta[2]<-Destimates[[1]]["m_Delta"]
  tau_Delta[2]<-1/Destimates[[1]]["sigma2_Delta"]
  alpha_CLINIC[2]<-Destimates[[2]]["alpha",1]
  beta_CLINIC[2]<-Destimates[[2]]["beta",1]
  alpha_HOME[2]<-Destimates[[2]]["alpha",2]
  beta_HOME[2]<-Destimates[[2]]["beta",2]
  
  rm(cleaned_nhanes_list,Sestimates,Destimates)
  
  xhat_i<- array(0, dim = c(4,2))
  xhat_f<- array(1, dim = c(4,2))
  resultslist<-list()
  
  # Convergence iterative looping over xhat values (where xhat is used to center all individuals by a 'mean'-ish value)
  while (sum(abs(xhat_i-xhat_f)/xhat_f)>0.005){
    print(paste0("Starting iteration ",toString(iii)))
    # Check for outputs of previous iterations
    filez<-paste0(savelocation, "exp_betas", run_name_string, mortality, ifelse(is.null(FR_score), "", paste0("_FRS",FR_score)) ,toString(iii-1),'.Rdata')
    print(" ")
    print(filez)
    print(" ")
    if(file.exists(filez)){
      print("Found the file exp_betas:")
      print(filez)
      file<-try(load(filez),silent=T)
      if(class(file)=="try-error") resultslist<-readRDS(file = filez)
      # Extract centering values from previous iteration
      M_av <- as.double(resultslist$xhat1)
      Delta_av <- as.double(resultslist$xhat2)
      tauis_av <- c(as.double(resultslist$xhat3),as.double(resultslist$xhat4))

      if(is.null(FR_score)){
        ixhat<-unname(rbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
      } else {
        ixhat<-unname(rbind(c(FRSc,0),resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
      }
      
      sumtable <- resultslist$fit.summary
      
      # Extract mean BP parameter values from previous iteration
      M_i<-cbind(sumtable[grep("M_i.*,1]",row.names(sumtable)),1],
                 sumtable[grep("M_i.*,2]",row.names(sumtable)),1])
      Delta_i<-cbind(sumtable[grep("Delta_i.*,1]",row.names(sumtable)),1],
                     sumtable[grep("Delta_i.*,2]",row.names(sumtable)),1])
      tauis_i<-cbind(sumtable[grep("tau.*,1]",row.names(sumtable)),1],
                     sumtable[grep("tau.*,2]",row.names(sumtable)),1],
                     sumtable[grep("tau.*,3]",row.names(sumtable)),1],
                     sumtable[grep("tau.*,4]",row.names(sumtable)),1])      
      
      rm(sumtable,resultslist)
    }  else{
      # If no previous iteration is found, generate the parameters from the BP measurements
      # Standard deviations of individuals 
      sds1 <- max(apply(xS, 1, sd)) # Systolic BP at clinic
      sds2 <- max(apply(xhomeS, 1, sd)) # Systolic BP at patients home
      sds3 <- max(apply(xD, 1, sd)) # Diastolic BP at clinic
      sds4 <- max(apply(xhomeD, 1, sd)) # Diastolic BP at patients home
      # Convert to precision
      tauis_i <- cbind(rep(1/sds1^2, N) , rep(1/sds2^2, N), rep(1/sds3^2, N), rep(1/sds4^2, N))
      # Calculate xhat centering values for precision
      tauis_av <- c(mean(1/sds1^2) , mean(1/sds2^2), mean(1/sds3^2), mean(1/sds4^2))
      # 'mean' BP between clinic and home
      M_i <- t(rbind(rowMeans( as.matrix(xS) ),rowMeans( as.matrix(xD) )))
      # Difference between clinic and home
      Delta_i <- t(rbind(rowMeans( abs((as.matrix(xS)-as.matrix(xhomeS))/2 ) ),rowMeans( abs((as.matrix(xD)-as.matrix(xhomeD))/2 ) )))
      # Calculate xhat centering values clinic and home mean and difference
      M_av <- c(mean(rowMeans( as.matrix(xS) )),mean(rowMeans( as.matrix(xD) )))
      Delta_av <- c(mean( rowMeans( abs((as.matrix(xS)-as.matrix(xhomeS))/2 ) )),mean( rowMeans( abs((as.matrix(xD)-as.matrix(xhomeD))/2 ) )))

      if(is.null(FR_score)){
        ixhat<-rbind(M_av,Delta_av,tauis_av[1:2],tauis_av[3:4])
      } else {
        ixhat<-rbind(c(FRSc,0),Delta_av,tauis_av[1:2],tauis_av[3:4])
      }
        
    }     
    
    if(!is.null(FR_score)) M_av<-c(FRSc,0)
    
    # Create lists of individual data to be read by rstan code
    standata <- list(
      Surv_T_and_delta = Surv_T_and_delta,
      N = N,
      max_nj = max(nj),
      yhomeS = as.matrix(xhomeS),
      yclinicS = as.matrix(xS) ,
      yhomeD = as.matrix(xhomeD),
      yclinicD = as.matrix(xD) ,    
      sex=sex,
      race=race,
      age= age,
      
      m_M=m_M,
      tau_M=tau_M,
      m_Delta=m_Delta,
      tau_Delta=tau_Delta,
      
      alpha_CLINIC=alpha_CLINIC,
      beta_CLINIC=beta_CLINIC,
      alpha_HOME=alpha_HOME,
      beta_HOME=beta_HOME,
      
      xhat_M=M_av,
      xhat_D=Delta_av,
      xhat_t=tauis_av  
      
    )
    # If using FRS value, add it to the list
    if(!is.null(FR_score)){standata<-c(standata, FRS=FRS)}
    
    # Initial guesses of parameters for MCMC 
    # my_inits <- function(chainnum){
    #   list(
    #     M_i = as.matrix(M_i),
    #     Delta_i = as.matrix(Delta_i),
    #     tauis = as.matrix(tauis_i)
    #   ) }
    my_inits <- function(chain_id){
      list(
        M_i = as.matrix(M_i),
        Delta_i = as.matrix(Delta_i),
        tauis = as.matrix(tauis_i)
      ) }
    
    # Print out to check the xhat centering values correspond to the data... FML    
    print("")
    print(c(mean(M_i[,1]),mean(M_i[,2])))
    print(c(mean(Delta_i[,1]),mean(Delta_i[,2])))
    print(c(mean(tauis_i[,1]),mean(tauis_i[,2]),mean(tauis_i[,3]),mean(tauis_i[,4])))
    print("")
    print(standata$alpha_CLINIC)
    print(standata$alpha_HOME)
    print(standata$beta_CLINIC)
    print(standata$beta_HOME)
    print("")
    print(mean(standata$yclinicS))
    print(mean(standata$yhomeS))
    print(mean(standata$yclinicD))
    print(mean(standata$yhomeD))
    print("")      
    print(standata$m_M)
    print(standata$tau_M)
    print(standata$m_Delta)
    print(standata$tau_Delta)
    print("")
    print(standata$xhat_M)  
    print(standata$xhat_D)  
    print(standata$xhat_t)        
    print("")
    ##################################################
    
    filename_to_save <- paste0(savelocation, "savedRun_", run_name_string, mortality,
                               ifelse(is.null(FR_score), "", paste0("_FRS",FR_score)),
                               "_",toString(iii),".csv")
    
    # Compile rstan code    
    # my_compiled_stan <- stan_model(file=my_file_stan)
    model<-cmdstan_model("./Model/mystanmodelFRS_DS_save.stan",compile = TRUE)
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN STAN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    # hwpstan<- sampling(my_compiled_stan, data=standata, init=my_inits , iter=50,
    #                    warmup=40, cores=mcs, chains=mcs, sample_file=filename_to_save)
    hwpstan<- model$sample(data=standata, init=my_inits , iter_warmup = 40,
                          iter_sampling = 50, chains = 1, parallel_chains = 1)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    
    print(paste0("Iteration ", iii," finito"))
    
    ############################## EXTRACT XHAT VALUES ######################################
    xhat1 <- xhat2 <- xhat3 <- xhat4 <- array(dim=c(mcs,2))
    NLL<-array(dim=mcs)
    
    # Loop over the number of MCMC chains used
    for(j in 1:mcs ){
      filename_to_save <- paste0(savelocation, "savedRun_", run_name_string, mortality,
                                 ifelse(is.null(FR_score), "", paste0("_FRS",FR_score)),
                                 "_",toString(iii),"_",toString(j),".csv")
      fit <- read_stan_csv(filename_to_save)
      
      ######################## Extract/generate parameters from chain ###########################
      beta<-as.matrix(fit,pars=c("beta"))
      gompertz<-as.matrix(fit,pars=c("gompertz"))
      gompz<-array(gompertz, c(dim(gompertz)[1],2,3,2))
      
      B<-theta<-array(0,dim=c(length(beta[,1]),N))
      for (k in 1:N){
        B[,k]<-gompz[,sex[k],race[k],1]
        theta[,k]<-gompz[,sex[k],race[k],2]
      }
      
      MMM<-as.matrix(fit,pars=c("M_i"))
      Delt<-as.matrix(fit,pars=c("Delta_i"))
      tauis<-as.matrix(fit,pars=c("tauis"))
      # Calculate summary values of parameters for next iteration
      fit.summary <- summary(fit)
      fit.summary <- fit.summary$summary
      # Calculate the Bayes factor
      BFs <- bayesFac(extract(fit, pars=c("beta")))
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
      
      xhat1[j,] <- xhat2[j,] <- xhat3[j,] <- xhat4[j,] <- 0.
      xhat<-ixhat
      
      print(ixhat)
      print(dim(ixhat))
      
      # xhat is defined differently based on FRS or nhanesA(inc. FRSpop)
      if(is.null(FR_score)){
        
        # xhat<-array(c(128,5.3,0.075,0.115,  75,3.85,0.15,0.19),dim = c(4,2))
        # xhat <-array(c(mean(M_i_S),mean(D_i_S),mean(tau_C_S),mean(tau_H_S),
        #               mean(M_i_D),mean(D_i_D),mean(tau_C_D),mean(tau_H_D)),dim = c(4,2))
          
        print(" ")
        print("New beta values")
        print("Systolic:")
        print(mean(beta[,1]))
        print(mean(beta[,2]))
        print(mean(beta[,3]))
        print(mean(beta[,4]))
        print("Diastolic:")
        print(mean(beta[,5]))
        print(mean(beta[,6]))
        print(mean(beta[,7]))
        print(mean(beta[,8]))    
        print(" ")
        
        ################### Newton Raphson method for 1D parameter root-solving #######################
        # I prefer to be explicit with defining the function and derivative so as not to make a mistake... FML
        
        func <- function(xhat) Survival_NLL_NF(M_i_S=M_i_S,M_i_D=M_i_D,D_i_S=D_i_S,D_i_D=D_i_D,
                                             tau_C_S=tau_C_S,tau_C_D=tau_C_D,tau_H_S=tau_H_S,tau_H_D=tau_H_D,
                                             age=age,delta=thisevent,T=T,beta=beta,B=B,theta=theta,xhat=xhat)
        
        output<-optim(ixhat,func)
        xhat<-output$par
        xhat[2,]%<>%abs
        
      }
      else{
        
        # xhat<-array(c(0,5.3,0.075,0.115,  0,3.85,0.15,0.19),dim = c(4,2))
        # xhat<-array(c(0,mean(D_i_S),mean(tau_C_S),mean(tau_H_S),
        #               0,mean(D_i_D),mean(tau_C_D),mean(tau_H_D)),dim = c(4,2))
          
        print(" ")
        print("New FRS beta values")
        print("FRS:")
        print(mean(beta[,1]))
        print("Systolic:")
        print(mean(beta[,2]))
        print(mean(beta[,3]))
        print(mean(beta[,4]))
        print("Diastolic:")
        print(mean(beta[,5]))
        print(mean(beta[,6]))
        print(mean(beta[,7]))
        print(" ")
        
        func <- function(xhat) Survival_NLL_F(FRS=FRS,FRSc=FRSc,D_i_S=D_i_S,D_i_D=D_i_D,
                                             tau_C_S=tau_C_S,tau_C_D=tau_C_D,tau_H_S=tau_H_S,tau_H_D=tau_H_D,
                                             age=age,delta=thisevent,T=T,beta=beta,B=B,theta=theta,xhat=xhat)
        
        output<-optim(ixhat[-1,],func)
        xhat<-ixhat
        xhat[-1,]<-output$par
        xhat[2,]%<>%abs
        # xhat[1,1]<-FRSc
        # xhat[1,2]<-0
        
      } 
      
      xhat1[j,]<-xhat[1,]
      xhat2[j,]<-xhat[2,]
      xhat3[j,]<-xhat[3,]
      xhat4[j,]<-xhat[4,]
      
      print(paste0("xhat values for MCMC chain: ",j))
      print(xhat1[j,]) # Should be zero for FRS
      print(xhat2[j,])
      print(xhat3[j,])
      print(xhat4[j,]) 

      NLL[j]<-output$value
        
      print(paste0("NLL[j]: ",NLL[j]))
      print(" ")
        
    }  

    resultslist <- list()

    ind<-which.min(NLL)
    NLL_n<-min(NLL)
    
    resultslist$xhat1<-xhat1[ind,]
    resultslist$xhat2<-xhat2[ind,]
    resultslist$xhat3<-xhat3[ind,]
    resultslist$xhat4<-xhat4[ind,]
    
    resultslist$NLL<-NLL_n
    
    # resultslist$xhat1<-colSums(xhat1*NLL/(sum(NLL)))
    # resultslist$xhat2<-colSums(xhat2*NLL/(sum(NLL)))
    # resultslist$xhat3<-colSums(xhat3*NLL/(sum(NLL)))
    # resultslist$xhat4<-colSums(xhat4*NLL/(sum(NLL)))
    
    print("xhat values:")
    print(resultslist$xhat1) # Should be zero for FRS
    print(resultslist$xhat2)
    print(resultslist$xhat3)
    print(resultslist$xhat4)
    print(" ")

    filename_to_save <- paste0(savelocation, "savedRun_", run_name_string, mortality,
                               ifelse(is.null(FR_score), "", paste0("_FRS",FR_score)),
                               "_",toString(iii),"_",toString(ind),".csv")
    fit <- read_stan_csv(filename_to_save)
      
    ######################## Extract/generate parameters from chain ###########################
    beta<-as.matrix(fit,pars=c("beta"))
    gompertz<-as.matrix(fit,pars=c("gompertz"))
    MMM<-as.matrix(fit,pars=c("M_i"))
    Delt<-as.matrix(fit,pars=c("Delta_i"))
    tauis<-as.matrix(fit,pars=c("tauis"))
    # Calculate summary values of parameters for next iteration
    fit.summary <- summary(fit)
    fit.summary <- fit.summary$summary
    # Calculate the Bayes factor
    BFs <- bayesFac(extract(fit, pars=c("beta")))
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
    
    resultslist$beta<-beta
    resultslist$gompertz<-gompertz      
    resultslist$BF<-BFs
    resultslist$fit.summary <-fit.summary

    resultslist$M_i_S<-M_i_S
    resultslist$M_i_D<-M_i_D
    resultslist$D_i_S<-D_i_S
    resultslist$D_i_D<-D_i_D
    resultslist$tau_H_S<-tau_H_S
    resultslist$tau_C_S<-tau_C_S
    resultslist$tau_H_D<-tau_H_D
    resultslist$tau_C_D<-tau_C_D
    
    save(resultslist, file=paste0(savelocation, "exp_betas", run_name_string, mortality, ifelse(is.null(FR_score), "", paste0("_FRS",FR_score)) ,toString(iii),'.Rdata'))
    
    xhat_i <- xhat_f
    xhat_f <- rbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4)
    
    # Check convergence parameter - don't vary more than 5%
    if(sum(abs(xhat_i-xhat_f)/xhat_f)>0.05) {rm(hwpstan)}
    
    iii<- iii+1
  }
  return(hwpstan)
}

############################ Load BP data #################################
# load(paste0(directory,"Data_cleaned/nhanes_cleaned_frames.RData"))
load(paste0(directory,"Data_cleaned/nhanes_cleaned_lists.RData"))

# sex=1-nhanes_passed$female
# Sex:
# 0 = female
# 1 = male

# race=nhanes_passed$white+2*nhanes_passed$mexican
# race:
# 0 = black
# 1 = white
# 2 = mexican

# type=sex*3+race+1
# To decode the type
# Remember that in DW's code, we removed "other" and mexican and other is the same
# (female)*3 + (black) + 1 
# (   0  )*3 + (  0  ) + 1 = 1
# (female)*3 + (white) + 1 
# (   0  )*3 + (  1  ) + 1 = 2
# (female)*3 + (mexic) + 1 
# (   0  )*3 + (  2  ) + 1 = 3

# ( male )*3 + (black) + 1 
# (   1  )*3 + (  0  ) + 1 = 4
# ( male )*3 + (white) + 1 
# (   1  )*3 + (  1  ) + 1 = 5
# ( male )*3 + (mexic) + 1 
# (   1  )*3 + (  2  ) + 1 = 6

# `mortality` can be:
#    "eventall"
#    "eventCVDHrt"
#    used to be one of these c("eventhrt", "eventother","eventCVD") but we don't have gompertz for them yet
# `FR_score` can be:
#     "ATP"
#     "1998" 
load(paste0(directory,"Model/BP_parameters.RData")) # Parameter estimates now BP_parameters
###################### Sort out the BP_parameters! #########################
# Function to put the parameters in the correct format for MCMC calculation
gamma_dimnames <- list( c('alpha','theta','beta') , c('Clinic' , 'Home'))
norm_dimnames <- c('m_M','m_Delta', 'sigma2_M', 'sigma2_Delta')
convert_param1 <- function(all_params){
  norm_params <- lapply(all_params, function(params) {
    with(params$Normal, c(m_M = M$m , m_Delta = Delta$m, sigma2_M = M$sigma2, sigma2_Delta = Delta$sigma2))
  })
  gamma_params <- lapply(all_params, function(params) {
    GP <- array( 0, dim = c(3,2) , dimnames = gamma_dimnames)
    GP[ 'alpha', ] <- c(params$Gamma$Clinic$alpha, params$Gamma$Home$alpha)
    GP[ 'theta', ] <- c(params$Gamma$Clinic$theta, params$Gamma$Home$theta)
    GP[ 'beta', ] <- GP[ 'alpha', ] / GP[ 'theta', ]
    GP
  })
  list(Normal = norm_params, Gamma = gamma_params )
}
estimates<-convert_param1(BP_parameters)
Sestimates=list(estimates$Normal$Systolic,estimates$Gamma$Systolic)    # Systolic BP estimates
Destimates=list(estimates$Normal$Diastolic,estimates$Gamma$Diastolic)  # Diastolic BP estimates

# rstan code
stan_file_no_frs <- paste0(directory,"Model/mystanmodel_DS.stan")
stan_file_with_frs <- paste0(directory,"Model/mystanmodelFRS_DS.stan")
# Which mortality to use - cardio-vascular and heart disease or all?
mortality <- c("eventCVDHrt", "eventall") 
# Which dataset?
nhaneslists <- c( "list_nhanesA", "list_nhanesFRS")

# Make a list to control all the parallel executions
par_control_list <- list(
list(runname="nhanesA_", nhaneslist= list_nhanesA,   mortality= "eventCVDHrt", stanfile=stan_file_no_frs, FR_score=NULL),
list(runname="nhanesA_", nhaneslist= list_nhanesA,   mortality=    "eventall", stanfile=stan_file_no_frs, FR_score=NULL),

list(runname="nhanesFRSpop_", nhaneslist= list_nhanesFRS, mortality= "eventCVDHrt", stanfile=stan_file_no_frs, FR_score=NULL),
list(runname="nhanesFRSpop_", nhaneslist= list_nhanesFRS, mortality=    "eventall", stanfile=stan_file_no_frs, FR_score=NULL),

list(runname="nhanesFRS_", nhaneslist= list_nhanesFRS, mortality= "eventCVDHrt", stanfile=stan_file_with_frs, FR_score="ATP"),
list(runname="nhanesFRS_", nhaneslist= list_nhanesFRS, mortality=    "eventall", stanfile=stan_file_with_frs, FR_score="ATP"),

list(runname="nhanesFRS_", nhaneslist= list_nhanesFRS, mortality= "eventCVDHrt", stanfile=stan_file_with_frs, FR_score="1998"),
list(runname="nhanesFRS_", nhaneslist= list_nhanesFRS, mortality=    "eventall", stanfile=stan_file_with_frs, FR_score="1998")
)

# Select the right data
pair <- par_control_list[[runnum]]

hwpstan<- RunNhanesModel(pair$stanfile,
                pair$nhaneslist,
                Sestimates,
                Destimates,
                pair$runname,
                pair$mortality,
                FR_score=pair$FR_score,
                iii=iii,
                mcs=mcs)

#save(hwpstan, file = paste0(directory,"Samples/HWPSTAN_DS_",toString(runnum),".Rdata"))

# Graph plotting and convergence checking
#library(shinystan)
#launch_shinystan(hwpstan)
