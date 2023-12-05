###########################################################################
##################### Linear Predictor Calculation ########################
###########################################################################

# Calculate the centering parameter automatically from the MCMC samples for the linear predictor term of the Gompertz equation
fxhat<-function(x,beta) {
  if(all(beta==0)) return(x)
  if(length(dim(x))>1) return((log(rowSums(exp(x*beta))) - log(ncol(x)) ) / beta)
  else return((log(sum(exp(x*beta))) - log(length(x)) ) / beta)
}

# The survival model which is used to predict the mortality risk of the individual.
# The term 'NF' refers to not-FRS whereby the mean diastolic and systolic blood pressure 
# are used instead of the FRS score of the individuals. 'F' is when the FRS score is used.
Survival_NLL_NF<-function(M_i_S,M_i_D,D_i_S,D_i_D,tau_C_S,tau_C_D,tau_H_S,tau_H_D,
                          age,delta,Time,beta,B,theta,xhat=NULL,roc=F,cumH=T){
  
  dimmy<-dim(D_i_S)
  
  # beta[ 1:1000 ]
  # x[ 1:1000  , 1:18018 ]
  if(is.null(xhat)){
    A1<-beta[,1]*(M_i_S       - fxhat(M_i_S,beta[,1]))
    A2<-beta[,2]*(abs(D_i_S)  - fxhat(abs(D_i_S),beta[,2]))
    A3<-beta[,3]*(tau_C_S     - fxhat(tau_C_S,beta[,3]))
    A4<-beta[,4]*(tau_H_S     - fxhat(tau_H_S,beta[,4]))
    A5<-beta[,5]*(M_i_D       - fxhat(M_i_D,beta[,5]))
    A6<-beta[,6]*(abs(D_i_D)  - fxhat(abs(D_i_D),beta[,6]))
    A7<-beta[,7]*(tau_C_D     - fxhat(tau_C_D,beta[,7]))
    A8<-beta[,8]*(tau_H_D     - fxhat(tau_H_D,beta[,8]))  
    
    linpred<- A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8
    
    if(roc){
      if(cumH) {H_j<-colMeans((B/theta)*exp(age*theta + linpred)*(exp(theta*Time) - 1))
      } else H_j<-colMeans(B*exp(theta*(age + Time) + linpred))
      d_j<-delta
      return(list(survy=data.frame(H_j=H_j,d_j=d_j, linpred=colMeans(linpred))))
    }
    
    n<-30
    breaks<-seq(from=min(age,na.rm = T)-1, to=max(age,na.rm = T)+1,length.out=n)
    
    a_j<-d_j<-vector(length = (n-2))
    H_j<-matrix(nrow = nrow(beta),ncol = (n-2))
    
    for (i in 2:(n-1)){
      
      # Age histogram a_j
      a_j[i-1]<-breaks[i]
      # indices of people that lived through age a_j
      ind<-age<=breaks[i]
      
      # cumulative hazard of all people that lived through age a_j
      if(cumH) {
        H_j[,i-1]<-rowSums((B[,ind]/theta[,ind])*exp(age[ind]*theta[,ind] + 
                                                       linpred[,ind])*(exp(theta[,ind]*Time[ind]) - 1),na.rm = T)
      } else {
        H_j[,i-1]<-rowSums((B[,ind])*exp(theta[,ind]*(age[ind] + Time[ind]) + linpred[,ind]),na.rm = T)
      }
      # cumulative deaths of all people that lived up to age a_j
      d_j[i-1]<-sum(delta[ind])
    }
    
    # Output the MLE estimate of xhat
    # (minimise the end as well as the mean distance between cumulative deaths - hazard)
    lenny<-ncol(H_j)
    iwin<-which.min(abs(H_j[,lenny]-d_j[lenny])*rowMeans(abs(H_j-d_j)))
    win<-abs(H_j[iwin,lenny]-d_j[lenny])*mean(abs(H_j[iwin,]-d_j))
    
    xhat<-array(dim=c(4,2))
    
    xhat[1,1]<-fxhat(M_i_S[iwin,],beta[iwin,1])
    xhat[2,1]<-fxhat(abs(D_i_S[iwin,]),beta[iwin,2])
    xhat[3,1]<-fxhat(tau_C_S[iwin,],beta[iwin,3])
    xhat[4,1]<-fxhat(tau_H_S[iwin,],beta[iwin,4])
    
    xhat[1,2]<-fxhat(M_i_D[iwin,],beta[iwin,5])
    xhat[2,2]<-fxhat(abs(D_i_D[iwin,]),beta[iwin,6])
    xhat[3,2]<-fxhat(tau_C_D[iwin,],beta[iwin,7])
    xhat[4,2]<-fxhat(tau_H_D[iwin,],beta[iwin,8])
    
    return(list(survy=data.frame(H_j=H_j[iwin,],d_j=d_j,a_j=a_j),xhat=xhat,iwin=iwin,win=win))
    
  } else {
    
    # xhat%<>%exp()%>%array(dim=c(4,2))
    
    A1<-beta[,1]*(M_i_S - xhat[1,1])
    A2<-beta[,2]*(abs(D_i_S) - xhat[2,1])
    A3<-beta[,3]*(tau_C_S - xhat[3,1])
    A4<-beta[,4]*(tau_H_S - xhat[4,1])
    A5<-beta[,5]*(M_i_D - xhat[1,2])
    A6<-beta[,6]*(abs(D_i_D) - xhat[2,2])
    A7<-beta[,7]*(tau_C_D - xhat[3,2])
    A8<-beta[,8]*(tau_H_D - xhat[4,2])  
  }
  linpred<- A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8
  
  if(roc){
    if(cumH) {H_j<-colMeans((B/theta)*exp(age*theta + linpred)*(exp(theta*Time) - 1))
    } else H_j<-colMeans(B*exp(theta*(age + Time) + linpred))
    d_j<-delta
    return(list(survy=data.frame(H_j=H_j,d_j=d_j, linpred=colMeans(linpred))))
  }
  
  n<-30
  breaks<-seq(from=min(age,na.rm = T)-1, to=max(age,na.rm = T)+1,length.out=n)
  
  a_j<-d_j<-vector(length = (n-2))
  H_j<-matrix(nrow = nrow(beta),ncol = (n-2))
  
  for (i in 2:(n-1)){
    
    # Age histogram a_j
    a_j[i-1]<-breaks[i]
    # indices of people that lived through age a_j
    ind<-age<=breaks[i]
    # cumulative hazard of all people that lived through age a_j
    if(cumH) {
      H_j[,i-1]<-rowSums((B[,ind]/theta[,ind])*exp(age[ind]*theta[,ind] + 
                                                     linpred[,ind])*(exp(theta[,ind]*Time[ind]) - 1),na.rm = T)
    } else {
      H_j[,i-1]<-rowSums((B[,ind])*exp(theta[,ind]*(age[ind] + Time[ind]) + linpred[,ind]),na.rm = T)
    }
    # cumulative deaths of all people that lived up to age a_j
    d_j[i-1]<-sum(delta[ind])
    
  }
  
  # Output the MLE estimate of xhat
  # (minimise the end as well as the mean distance between cumulative deaths - hazard)
  lenny<-ncol(H_j)
  iwin<-which.min(abs(H_j[,lenny]-d_j[lenny])*rowMeans(abs(H_j-d_j)))
  win<-abs(H_j[iwin,lenny]-d_j[lenny])*mean(abs(H_j[iwin,]-d_j))
  
  return(list(survy=data.frame(H_j=sort(H_j[iwin,]),d_j=d_j,a_j=a_j),iwin=iwin,win=win))
  
}

# Using the FRS score instead of the mean systolic & diastolic blood pressure
Survival_NLL_F<-function(FRS,FRSc,D_i_S,D_i_D,tau_C_S,tau_C_D,tau_H_S,tau_H_D,
                         age,delta,Time,beta,B,theta,xhat=NULL,roc=F,cumH=T){
  
  dimmy<-dim(D_i_S)
  
  # beta[ 1:1000 ]
  # x[ 1:1000  , 1:18018 ]
  if(is.null(xhat)){
    A1<-beta[,1]%o%(FRS-FRSc)
    A2<-beta[,2]*(abs(D_i_S)  - fxhat(abs(D_i_S),beta[,2]))
    A3<-beta[,3]*(tau_C_S     - fxhat(tau_C_S,beta[,3]))
    A4<-beta[,4]*(tau_H_S     - fxhat(tau_H_S,beta[,4]))
    A5<-beta[,5]*(abs(D_i_D)  - fxhat(abs(D_i_D),beta[,5]))
    A6<-beta[,6]*(tau_C_D     - fxhat(tau_C_D,beta[,6]))
    A7<-beta[,7]*(tau_H_D     - fxhat(tau_H_D,beta[,7]))  
    
    linpred<- A1 + A2 + A3 + A4 + A5 + A6 + A7
    
    if(roc){
      if(cumH) {H_j<-colMeans((B/theta)*exp(age*theta + linpred)*(exp(theta*Time) - 1))
      } else H_j<-colMeans(B*exp(theta*(age + Time) + linpred))
      d_j<-delta
      return(list(survy=data.frame(H_j=H_j,d_j=d_j, linpred=colMeans(linpred))))
    }
    
    n<-30
    breaks<-seq(from=min(age,na.rm = T)-1, to=max(age,na.rm = T)+1,length.out=n)
    
    a_j<-d_j<-vector(length = (n-2))
    H_j<-matrix(nrow = nrow(beta),ncol = (n-2))
    
    for (i in 2:(n-1)){
      
      # Age histogram a_j
      a_j[i-1]<-breaks[i]
      # indices of people that lived through age a_j
      ind<-age<=breaks[i]
      
      # cumulative hazard of all people that lived through age a_j
      if(cumH) {
        H_j[,i-1]<-rowSums((B[,ind]/theta[,ind])*exp(age[ind]*theta[,ind] + 
                                                       linpred[,ind])*(exp(theta[,ind]*Time[ind]) - 1),na.rm = T)
      } else {
        H_j[,i-1]<-rowSums((B[,ind])*exp(theta[,ind]*(age[ind] + Time[ind]) + linpred[,ind]),na.rm = T)
      }
      # cumulative deaths of all people that lived up to age a_j
      d_j[i-1]<-sum(delta[ind])
    }
    
    # Output the MLE estimate of xhat
    # (minimise the end as well as the mean distance between cumulative deaths - hazard)
    lenny<-ncol(H_j)
    iwin<-which.min(abs(H_j[,lenny]-d_j[lenny])*rowMeans(abs(H_j-d_j)))
    win<-abs(H_j[iwin,lenny]-d_j[lenny])*mean(abs(H_j[iwin,]-d_j))
    
    xhat<-array(dim=c(4,2))
    
    xhat[1,1:2]<-c(FRSc,0)
    
    xhat[2,1]<-fxhat(abs(D_i_S[iwin,]),beta[iwin,2])
    xhat[3,1]<-fxhat(tau_C_S[iwin,],beta[iwin,3])
    xhat[4,1]<-fxhat(tau_H_S[iwin,],beta[iwin,4])
    
    xhat[2,2]<-fxhat(abs(D_i_D[iwin,]),beta[iwin,5])
    xhat[3,2]<-fxhat(tau_C_D[iwin,],beta[iwin,6])
    xhat[4,2]<-fxhat(tau_H_D[iwin,],beta[iwin,7])
    
    return(list(survy=data.frame(H_j=H_j[iwin,],d_j=d_j,a_j=a_j),xhat=xhat,iwin=iwin,win=win))
    
  } else {
    
    A1<-beta[,1]%o%(FRS-FRSc)
    A2<-beta[,2]*(abs(D_i_S) - xhat[2,1])
    A3<-beta[,3]*(tau_C_S - xhat[3,1])
    A4<-beta[,4]*(tau_H_S - xhat[4,1])
    A5<-beta[,5]*(abs(D_i_D) - xhat[2,2])
    A6<-beta[,6]*(tau_C_D - xhat[3,2])
    A7<-beta[,7]*(tau_H_D - xhat[4,2])
  }
  
  linpred<-A1 + A2 + A3 + A4 + A5 + A6 + A7
  
  if(roc){
    if(cumH) {H_j<-colMeans((B/theta)*exp(age*theta + linpred)*(exp(theta*Time) - 1))
    } else H_j<-colMeans(B*exp(theta*(age + Time) + linpred))
    d_j<-delta
    return(list(survy=data.frame(H_j=H_j,d_j=d_j, linpred=colMeans(linpred))))
  }
  
  n<-30
  breaks<-seq(from=min(age,na.rm = T)-1, to=max(age,na.rm = T)+1,length.out=n)
  
  a_j<-d_j<-vector(length = (n-2))
  H_j<-matrix(nrow = nrow(beta),ncol = (n-2))
  
  for (i in 2:(n-1)){
    
    # Age histogram a_j
    a_j[i-1]<-breaks[i]
    # indices of people that lived through age a_j
    ind<-age<=breaks[i]
    
    # cumulative hazard of all people that lived through age a_j
    if(cumH) {
      H_j[,i-1]<-rowSums((B[,ind]/theta[,ind])*exp(age[ind]*theta[,ind] + 
                                                     linpred[,ind])*(exp(theta[,ind]*Time[ind]) - 1),na.rm = T)
    } else {
      H_j[,i-1]<-rowSums((B[,ind])*exp(theta[,ind]*(age[ind] + Time[ind]) + linpred[,ind]),na.rm = T)
    }
    # cumulative deaths of all people that lived up to age a_j
    d_j[i-1]<-sum(delta[ind])
  }
  
  # Output the MLE estimate of xhat
  # (minimise the end as well as the mean distance between cumulative deaths - hazard)
  lenny<-ncol(H_j)
  iwin<-which.min(abs(H_j[,lenny]-d_j[lenny])*rowMeans(abs(H_j-d_j)))
  win<-abs(H_j[iwin,lenny]-d_j[lenny])*mean(abs(H_j[iwin,]-d_j))
  
  return(list(survy=data.frame(H_j=sort(H_j[iwin,]),d_j=d_j,a_j=a_j),iwin=iwin,win=win))
  
}

# Build a database for the Nhanes data
DF_Nhanes<-function(list_nhanes){
  
  DF<-as.data.frame.list(list_nhanes[c("eventall","eventCVDHrt","black","white","other",
                                       "female","male","T","age")])
  if(!is.null(list_nhanes[["FRS.ATP"]]) & !is.null(list_nhanes[["FRS.1998"]]))
    DF%<>%cbind(as.data.frame.list(list_nhanes[c("FRS.ATP","FRS.1998")]))
  
  return(DF)
}

FilterRL<-function(RL,ind,sigma=T){
  
  if(any(grepl(names(RL),pattern = "tau_"))){ columns<-c("M_i_S","M_i_D","D_i_S","D_i_D","tau_H_S","tau_C_S","tau_H_D","tau_C_D")
  } else {columns<-c("M_i_S","M_i_D","D_i_S","D_i_D","sigma_H_S","sigma_C_S","sigma_H_D","sigma_C_D")}
  
  for(c in columns){
    RL[[c]]<-RL[[c]][,ind]
  }
  
  return(RL)
  
}

ReduceDemog<-function(RL,Ethnicity=NULL,Gender=NULL,ageBnds=NULL,othDeaths=F, ind=NULL){
  
  if(RL$FRS) nhanes<-DF_Nhanes(list_nhanesFRS) else nhanes<-DF_Nhanes(list_nhanesA)
  
  if(is.null(ind) & !othDeaths) {ind<-rep(T,nrow(nhanes)) 
  }else if(is.null(ind) & othDeaths) {ind<-nhanes$eventall+nhanes$eventCVDHrt<2}
  
  if(!is.null(Ethnicity) | !is.null(Gender) | !is.null(ageBnds)){
    # First, ethnicity
    if(!is.null(Ethnicity)) ind<-ind & as.logical(nhanes[[Ethnicity]])
    if(!is.null(Gender)) ind<-ind & as.logical(nhanes[[Gender]])
    if(!is.null(ageBnds)) {
      # First check the bounds are correctly specified
      if(ageBnds[2]<ageBnds[1]) stop("misspecified age bounds to use in GetSurvival")
      ind<-ind & nhanes$age>=ageBnds[1] & nhanes$age<ageBnds[2]
    }
    nhanes<-nhanes[ind,]
    RL%<>%FilterRL(ind, sigma=T)
  }
  return(list(nhanes=nhanes,RL=RL))
}

# Make the predictions on the survival outcomes of the individuals, 
# based on posterior samples from the HMC algorithm (from Stan)
GetSurvival<-function(RL,roc=F,Ethnicity=NULL,Gender=NULL,ageBnds=NULL,Year=20,usexhat=TRUE, cumH=T, sigma=T,othDeaths=F, ind=NULL){
  
  if(usexhat) xhat<-t(cbind(RL$xhat1,RL$xhat2,RL$xhat3,RL$xhat4))
  
  if(sigma) variancer <- "sigma" else variancer <- "tau"
  
  if(RL$FRS) nhanes<-DF_Nhanes(list_nhanesFRS) else nhanes<-DF_Nhanes(list_nhanesA)
  
  tmp<-ReduceDemog(RL,Ethnicity,Gender,ageBnds,othDeaths = othDeaths, ind=ind); nhanes<-tmp$nhanes; RL<-tmp$RL; rm(tmp)
  # Set the year since start of study to calculate
  nhanes$T<-Year
  
  LLL<-length(nhanes$T)
  N<-dim(RL$beta)[1]
  B<-theta<-array(0,dim=c(N,LLL))
  
  gender<-nhanes$male+1
  ethn<-nhanes$black+2*nhanes$white+3*nhanes$other
  
  if(RL$eventall){delta<-nhanes$eventall
  } else{delta<-nhanes$eventCVDHrt}
  
  for (j in 1:LLL){
    B[,j]<-RL$gompertz[,gender[j],ethn[j],1]
    theta[,j]<-RL$gompertz[,gender[j],ethn[j],2]
  }
  rm(gender,ethn)
  
  if(!is.na(RL$FRSt)){
    
    if(RL$FRSt=='ATP'){FRS<-nhanes$FRS.ATP
    } else if(RL$FRSt=='1998'){FRS<-nhanes$FRS.1998
    } else{stop("Wrong FRS entry")}
    
    FRSc<-log(mean(exp(FRS)))
    
    if(usexhat)  {survy<-(Survival_NLL_F(FRS,FRSc,RL$D_i_S,RL$D_i_D,RL[[paste0(variancer,"_C_S")]],RL[[paste0(variancer,"_C_D")]],
                                         RL[[paste0(variancer,"_C_D")]],RL[[paste0(variancer,"_H_D")]],nhanes$age,delta,
                                         nhanes$T,RL$beta,B,theta,xhat,roc=roc,cumH=cumH))$survy
    } else {
      survy<-(Survival_NLL_F(FRS,FRSc,RL$D_i_S,RL$D_i_D,RL[[paste0(variancer,"_C_S")]],RL[[paste0(variancer,"_C_D")]],
                             RL[[paste0(variancer,"_C_D")]],RL[[paste0(variancer,"_H_D")]],nhanes$age,delta,
                             nhanes$T,RL$beta,B,theta,roc=roc,cumH=cumH))$survy
    }
    rm(FRS)
    
  } else{
    
    if(usexhat)  {survy<-(Survival_NLL_NF(RL$M_i_S,RL$M_i_D,RL$D_i_S,RL$D_i_D,RL[[paste0(variancer,"_C_S")]],RL[[paste0(variancer,"_C_D")]],
                                          RL[[paste0(variancer,"_C_D")]],RL[[paste0(variancer,"_H_D")]],nhanes$age,delta,
                                          nhanes$T,RL$beta,B,theta,xhat,roc=roc,cumH=cumH))$survy
    } else {
      survy<-(Survival_NLL_NF(RL$M_i_S,RL$M_i_D,RL$D_i_S,RL$D_i_D,RL[[paste0(variancer,"_C_S")]],RL[[paste0(variancer,"_C_D")]],
                              RL[[paste0(variancer,"_C_D")]],RL[[paste0(variancer,"_H_D")]],nhanes$age,delta,
                              nhanes$T,RL$beta,B,theta,roc=roc,cumH=cumH))$survy
    }
    
  }
  
  return(survy) 
  
}

# Normalise the beta values
ConvertBeta<-function(RL,FRS=NULL){
  
  N<-nrow(RL$beta)
  out<-data.frame()
  if(is.null(FRS)){
    out%<>%rbind(data.frame(variable=rep("M_i_S",N),value=RL$beta[,1],normalised=RL$beta[,1]*apply(RL$M_i_S,1,sd)))
    out%<>%rbind(data.frame(variable=rep("D_i_S",N),value=RL$beta[,2],normalised=RL$beta[,2]*apply(RL$D_i_S,1,sd)))
    out%<>%rbind(data.frame(variable=rep(paste0(variancer,"_C_S"),N),value=RL$beta[,3],normalised=RL$beta[,3]*apply(RL[[paste0(variancer,"_C_S")]],1,sd)))
    out%<>%rbind(data.frame(variable=rep(paste0(variancer,"_H_S"),N),value=RL$beta[,4],normalised=RL$beta[,4]*apply(RL[[paste0(variancer,"_H_S")]],1,sd)))
    
    out%<>%rbind(data.frame(variable=rep("M_i_D",N),value=RL$beta[,5],normalised=RL$beta[,5]*apply(RL$M_i_D,1,sd)))
    out%<>%rbind(data.frame(variable=rep("D_i_D",N),value=RL$beta[,6],normalised=RL$beta[,6]*apply(RL$D_i_D,1,sd)))
    out%<>%rbind(data.frame(variable=rep(paste0(variancer,"_C_D"),N),value=RL$beta[,7],normalised=RL$beta[,7]*apply(RL[[paste0(variancer,"_C_D")]],1,sd)))
    out%<>%rbind(data.frame(variable=rep(paste0(variancer,"_H_D"),N),value=RL$beta[,8],normalised=RL$beta[,8]*apply(RL[[paste0(variancer,"_H_D")]],1,sd)))
    
  } else {
    out%<>%rbind(data.frame(variable=rep("FRS",N),value=RL$beta[,1],normalised=RL$beta[,1]*sd(FRS)))
    out%<>%rbind(data.frame(variable=rep("D_i_S",N),value=RL$beta[,2],normalised=RL$beta[,2]*apply(RL$D_i_S,1,sd)))
    out%<>%rbind(data.frame(variable=rep(paste0(variancer,"_C_S"),N),value=RL$beta[,3],normalised=RL$beta[,3]*apply(RL[[paste0(variancer,"_C_S")]],1,sd)))
    out%<>%rbind(data.frame(variable=rep(paste0(variancer,"_H_S"),N),value=RL$beta[,4],normalised=RL$beta[,4]*apply(RL[[paste0(variancer,"_H_S")]],1,sd)))
    
    out%<>%rbind(data.frame(variable=rep("D_i_D",N),value=RL$beta[,5],normalised=RL$beta[,5]*apply(RL$D_i_D,1,sd)))
    out%<>%rbind(data.frame(variable=rep(paste0(variancer,"_C_D"),N),value=RL$beta[,6],normalised=RL$beta[,6]*apply(RL[[paste0(variancer,"_C_D")]],1,sd)))
    out%<>%rbind(data.frame(variable=rep(paste0(variancer,"_H_D"),N),value=RL$beta[,7],normalised=RL$beta[,7]*apply(RL[[paste0(variancer,"_H_D")]],1,sd)))
    
  }
  
  return(out)
  
}

# Normalise the gompertz parameters
ConvertGompertz<-function(gompertz){
  # Female <-gender[j]=2, male <-gender[j]=1
  # Black <- ethn[j]=1, white<-ethn[j]=2, other<-ethn[j]=3
  N<-dim(gompertz)[1]
  gompz<-data.frame()
  
  gompz%<>%rbind(data.frame(Gender=rep("Female",N),Ethnicity=rep("Black",N),
                            B=gompertz[,1,1,1],theta=gompertz[,1,1,2]))
  gompz%<>%rbind(data.frame(Gender=rep("Female",N),Ethnicity=rep("White",N),
                            B=gompertz[,1,2,1],theta=gompertz[,1,2,2]))
  gompz%<>%rbind(data.frame(Gender=rep("Female",N),Ethnicity=rep("Other",N),
                            B=gompertz[,1,3,1],theta=gompertz[,1,3,2]))
  gompz%<>%rbind(data.frame(Gender=rep("Male",N),Ethnicity=rep("Black",N),
                            B=gompertz[,2,1,1],theta=gompertz[,2,1,2]))
  gompz%<>%rbind(data.frame(Gender=rep("Male",N),Ethnicity=rep("White",N),
                            B=gompertz[,2,2,1],theta=gompertz[,2,2,2]))
  gompz%<>%rbind(data.frame(Gender=rep("Male",N),Ethnicity=rep("Other",N),
                            B=gompertz[,2,3,1],theta=gompertz[,2,3,2]))
  
  return(gompz)
  
}

### CALCULATE ROC FUNCTION ###
calculate_roc <- function(df, cost_of_fp=1, cost_of_fn=1, n=100) {
  tpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$survived == 1) / sum(df$survived == 1)
  }
  
  fpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$survived == 0) / sum(df$survived == 0)
  }
  
  cost <- function(df, threshold, cost_of_fp, cost_of_fn) {
    sum(df$pred >= threshold & df$survived == 0) * cost_of_fp + 
      sum(df$pred < threshold & df$survived == 1) * cost_of_fn
  }
  
  roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th))
  roc$cost <- sapply(roc$threshold, function(th) cost(df, th, cost_of_fp, cost_of_fn))
  
  dx <- diff(roc$fpr)
  my <- (roc$tpr[1:(n - 1)] + roc$tpr[2:n]) / 2
  roc$auroc<-sum(abs(dx*my))
  
  return(roc)
}
### PLOT THAT ROC FUNCTION ###
plot_roc <- function(roc, threshold, cost_of_fp, cost_of_fn, costy=TRUE,trad=T) {
  
  library(gridExtra)
  
  norm_vec <- function(v) (v - min(v))/diff(range(v))
  
  idx_threshold = which.min(abs(roc$threshold-threshold))
  
  col_ramp <- colorRampPalette(c("green","orange","red","black"))(100)
  col_by_cost <- col_ramp[ceiling(norm_vec(roc$cost)*99)+1]
  
  if(trad){
    p_roc <- ggplot(roc, aes(fpr,tpr)) + 
      geom_line(color=rgb(0,0,1,alpha=0.3)) +
      geom_point(color='red', size=4, alpha=0.5) +
      coord_fixed() + ylim(c(0,1)) +
      # geom_line(aes(threshold,threshold), color=rgb(0,0,1,alpha=0.5)) +
      labs(title = sprintf("ROC Curve")) + xlab("False Positive Rate") + ylab("True Positive Rate") + theme(plot.title = element_text(hjust = 0.5)) +
      #geom_hline(yintercept=roc[idx_threshold,"tpr"], alpha=0.5, linetype="dashed") +
      #geom_vline(xintercept=roc[idx_threshold,"fpr"], alpha=0.5, linetype="dashed") +
      annotate("text", x = .3, y = .1, size=5, label = paste0("AUROC =", round(roc$auroc, 2)))
  }else {  
    p_roc <- ggplot(roc, aes(fpr,(tpr-fpr))) + 
      geom_line(color=rgb(0,0,1,alpha=0.3),size=3) +
      geom_point(color=col_by_cost, size=4, alpha=0.5) +
      coord_fixed() + ylim(c(0,0.5)) +
      # geom_line(aes(threshold,threshold), color=rgb(0,0,1,alpha=0.5)) +
      labs(title = sprintf("ROC Curve")) + xlab("False Positive Rate") + ylab("True Positive Rate - False Positive Rate") + theme(plot.title = element_text(hjust = 0.5)) +
      #geom_hline(yintercept=roc[idx_threshold,"tpr"], alpha=0.5, linetype="dashed") +
      #geom_vline(xintercept=roc[idx_threshold,"fpr"], alpha=0.5, linetype="dashed") +
      annotate("text", x = .3, y = .1, size=5, label = paste0("AUROC =", round(roc$auroc, 2)))
  }
  
  if(costy){
    p_cost <- ggplot(roc, aes(threshold, cost)) +
      geom_line(color=rgb(0,0,1,alpha=0.3)) +
      geom_point(color=col_by_cost, size=4, alpha=0.5) +
      labs(title = sprintf("Cost Function")) + xlab("Threshold") + ylab("Cost") + theme(plot.title = element_text(hjust = 0.5))
    #geom_vline(xintercept=threshold, alpha=0.5, linetype="dashed")
    
    grid.arrange(p_roc, p_cost, ncol=2)
  } else {p_roc}
  #sub_title <- sprintf("threshold at %.2f - cost of FP = %d, cost of FN = %d", threshold, cost_of_fp, cost_of_fn)
  #, sub=textGrob(sub_title, gp=gpar(cex=1), just="bottom"))
}


# replace Time in survival_NLL functions
# Then spit out the FP & TPs, dependent upon the threshold value:
ConfVals<-function(df,Year,threshold){    
  if((df[["Time"]])<Year){
    if(df[["delta"]]==1) {
      # TP if model predicts alive and the censoring occurred after Year, else we cannot tell
      outcome<-ifelse((df[["pred"]]>=threshold),"TP","FN")
      # If there is no death then we cannot know whether the model correctly predicted or not
    } else outcome<-NA 
  } else if(df[["Time"]]>=Year){
    # FP if the model predicts death, TN otherwise
    outcome<-ifelse((df[["pred"]]>=threshold),"FP","TN")
  } else outcome<-NA 
  
  return(outcome)
}

calc_Year_roc<-function(RL,Year=20L,Ethnicity=NULL,Gender=NULL,ageBnds=NULL, RedCovars=NULL, othDeaths=F, ind=NULL){
  # Calculate the cumulative hazard values per individual at the given time since cohort start
  
  if(is.null(RedCovars)) {
    if(is.na(RL$FRSt)) RL$beta[,c(3,4,7,8)]<-0 else RL$beta[,c(3,4,6,7)]<-0
  } else if (RedCovars=="sysmean") {
    if(is.na(RL$FRSt)) RL$beta[,2:8]<-0 else RL$beta[,2:7]<-0
  } else if (RedCovars=="deltas") {
    if(is.na(RL$FRSt)) RL$beta[,c(1,3,4,5,7,8)]<-0 else RL$beta[,c(1,3,4,6,7)]<-0
  }
  
  survy<-GetSurvival(RL,roc=T,Ethnicity = Ethnicity,Gender = Gender,Year = Year,ageBnds = ageBnds, othDeaths=othDeaths, ind=ind)
  # Modify the nhanes and RL objects to reduce to the same demography
  tmp<-ReduceDemog(RL = RL,Ethnicity = Ethnicity,Gender = Gender,ageBnds = ageBnds,othDeaths = othDeaths, ind=ind); nhanes<-tmp$nhanes; RL<-tmp$RL; rm(tmp)
  eventer<-ifelse(RL$eventall,"eventall","eventCVDHrt")
  df<-data.frame(Time=nhanes$T,delta=nhanes[[eventer]],pred=survy$H_j)
  # How many values do you want to use to calculate the ROC & AUC values?
  n<-100
  
  roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
  PredPerf <- sapply(roc$threshold, function(th) apply(df,1,function(df) ConfVals(df,Year,th)))
  
  roc$tpr<-vapply(1:n,function(i) sum(PredPerf[,i]=="TP",na.rm = T)/sum((PredPerf[,i]=="FN" | PredPerf[,i]=="TP"),na.rm = T),
                  FUN.VALUE = numeric(1))
  roc$fpr<-vapply(1:n,function(i) sum(PredPerf[,i]=="FP",na.rm = T)/sum((PredPerf[,i]=="FP" | PredPerf[,i]=="TN"),na.rm = T),
                  FUN.VALUE = numeric(1))
  
  dx <- diff(roc$fpr)
  my <- (roc$tpr[1:(n - 1)] + roc$tpr[2:n]) / 2
  roc$auroc<-sum(abs(dx*my))
  
  return(roc)
}






















