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
# install.packages("rstan", repos="https://cloud.r-project.org",lib="../Rstan")
.libPaths()
library(pracma)
library(plyr)
library(magrittr)
library(rstan,lib.loc="../Rstan/")
rstan_options(auto_write = TRUE)

directory<- "/data/localhost/patten/BP_NhanesA/"

# Read the run number
args <- commandArgs(trailingOnly = TRUE)

runnum <- as.integer(args[1])
print(paste0("Calculating for run number: ",runnum))

iii<-0
mcs<-4
sigma <- FALSE
autopred <- FALSE
# Set the starting point for iterative simulations of the xhat parameter
if(length(args)>1) iii<- as.integer(args[2])
print(paste0("xhat iteration number: ",iii))

if(length(args)>2) {
    sigma<- as.logical(args[3])
}
if(sigma) print("Using sigma not tau for variance term") else print("Using tau for variance term")

if(length(args)>3) {
    autopred<- as.logical(args[4])
}
if(autopred) print("Automatically calculating xhat centering values") else print("Using fixed xhat centering values")

# Set the number of parallel computational tasks
if(length(args)>4) mcs<- as.integer(args[5])
print(paste0("Number of parallel cpu tasks: ",mcs))

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
                           sigma=FALSE,
                           autopred=FALSE,
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
  # run_name_string<-paste0(run_name_string,"_RD2")
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
  #while (sum(abs(xhat_i-xhat_f)/xhat_f)>0.005){
    print(paste0("Starting iteration ",toString(iii)))
    
    # Check for outputs of previous iterations
    ## filez<-paste0(savelocation, "exp_betas", run_name_string, mortality, ifelse(is.null(FR_score), "", paste0("_FRS",FR_score)) ,ifelse(sigma,"_sigma_","_tau_"),ifelse(autopred,"autopred_",""),toString(iii-1),'.Rdata')
    filez<-paste0(savelocation, "exp_betas", run_name_string, mortality, ifelse(is.null(FR_score), "", paste0("_FRS",FR_score)) ,ifelse(sigma,"_sigma_","_tau_"),toString(iii-1),'.Rdata')

    print(" ")
    print(filez)
    print(" ")
    
    if(file.exists(filez)){
      print("Found the file exp_betas:")
      print(filez)
      file<-try(load(filez),silent=T)
      if(class(file)=="try-error") resultslist<-readRDS(file = filez)
      # Extract centering values from previous iteration
      if(is.null(FR_score)){
        ixhat<-unname(rbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
      } else {
        ixhat<-unname(rbind(c(FRSc,0),resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
      }
      
      sumtable <- resultslist$fit.summary
      
      M_i <- cbind(colMeans(resultslist$M_i_S),colMeans(resultslist$M_i_D))
      Delta_i <- cbind(colMeans(resultslist$D_i_S),colMeans(resultslist$D_i_D))     

      if(sigma){
          tauis_i <- 1./sqrt(cbind(colMeans(resultslist$sigma_C_S),colMeans(resultslist$sigma_H_S),
                       colMeans(resultslist$sigma_C_D),colMeans(resultslist$sigma_H_D)))
      } else {
          tauis_i <- cbind(colMeans(resultslist$tau_C_S),colMeans(resultslist$tau_H_S),
                       colMeans(resultslist$tau_C_D),colMeans(resultslist$tau_H_D))
      }
      
      beta_i  <- colMeans(resultslist$beta)
      gompertz_i  <- colMeans(array(resultslist$gompertz,c(nrow(resultslist$gompertz),2,3,2)))

      my_inits <- function(chainnum){
          list(
              M_i = as.matrix(M_i),
              Delta_i = as.matrix(Delta_i),
              tauis = as.matrix(tauis_i),
              beta  =  as.vector(beta_i),
              gompertz = gompertz_i   
          )
      }
      
      rm(sumtable,resultslist)
    }
    else{
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
      sigmais_av <- 1/sqrt(tauis_av)  
      # 'mean' BP between clinic and home
      M_i <- t(rbind(rowMeans( as.matrix(xS) ),rowMeans( as.matrix(xD) )))
      # Difference between clinic and home
      Delta_i <- t(rbind(rowMeans( abs((as.matrix(xS)-as.matrix(xhomeS))/2 ) ),rowMeans( abs((as.matrix(xD)-as.matrix(xhomeD))/2 ) )))
      # Calculate xhat centering values clinic and home mean and difference
      M_av <- c(mean(rowMeans( as.matrix(xS) )),mean(rowMeans( as.matrix(xD) )))
        Delta_av <- c(mean( rowMeans( abs((as.matrix(xS)-as.matrix(xhomeS))/2 ) )),mean( rowMeans( abs((as.matrix(xD)-as.matrix(xhomeD))/2 ) )))

      my_inits <- function(chainnum){
          list(
              M_i = as.matrix(M_i),
              Delta_i = as.matrix(Delta_i),
              tauis = as.matrix(tauis_i)
          )
      }

      if(sigma){
          ixhat<-rbind(M_av,Delta_av,sigmais_av[1:2],sigmais_av[3:4])
      } else {
          ixhat<-rbind(M_av,Delta_av,tauis_av[1:2],tauis_av[3:4])
      }

    }

    if(!is.null(FR_score)) ixhat[1,]<-c(FRSc,0)
    
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
      beta_HOME=beta_HOME      
    )

    if(!autopred) standata<-c(standata,list(xhat_M=ixhat[1,],xhat_D=ixhat[2,],xhat_v=c(ixhat[3:4,])))        
    # If using FRS value, add it to the list
    if(!is.null(FR_score)) standata<-c(standata,list(FRSc=ixhat[1,]),FRS=FRS)    
    
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
    ##################################################
    
    filename_to_save <- paste0(savelocation, "savedRun_", run_name_string, mortality,
                               ifelse(sigma, "_sigma", "_tau"),
                               ifelse(is.null(FR_score), "", paste0("_FRS",FR_score)),
                               ifelse(autopred,"_autopred_","_"),
                               toString(iii),".csv")
    
    # Compile rstan code    
    my_compiled_stan <- stan_model(file=my_file_stan)
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN STAN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    hwpstan<- sampling(my_compiled_stan, data=standata, init=my_inits , iter=6000,
                       warmup=4500, cores=mcs, chains=mcs, sample_file=filename_to_save)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    
    print(paste0("Iteration ", iii," finito"))
    
}

############################ Load BP data #################################
load(paste0(directory,"Data_cleaned/nhanes_cleaned_frames.RData"))
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
# stan_file_no_frs <- paste0(directory,"Model/mystanmodel_DS_sigma_v2.stan")
# stan_file_with_frs <- paste0(directory,"Model/mystanmodelFRS_DS_sigma_v2.stan")

stan_file_no_frs <- paste0(directory,"Model/mystanmodel_DS_",ifelse(sigma,"sigma","tau"),"_v2",ifelse(autopred,"_autopred",""),".stan")
stan_file_with_frs <- paste0(directory,"Model/mystanmodelFRS_DS_",ifelse(sigma,"sigma","tau"),"_v2",ifelse(autopred,"_autopred",""),".stan")

# Which mortality to use - cardio-vascular and heart disease or all?
mortality <- c("eventCVDHrt", "eventall") 

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
                mcs=mcs,
                sigma=sigma,
                autopred=autopred)

save(hwpstan, file = paste0(directory,"Samples/HWPSTAN_DS_",toString(runnum),".Rdata"))
