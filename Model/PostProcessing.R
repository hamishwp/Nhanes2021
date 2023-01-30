library(tidyverse)
library(latex2exp)
library(ggplot2)
library(dplyr)
library(pracma)
library(magrittr)
library(ggpubr)
# devtools::install_github("ottlngr/LexisPlotR")
library(LexisPlotR)
source("./Model/Functions.R")

directory<-'~/Document/Coding/Oxford/Nhanes2021/'

saveresultslist<-function(resultslist,RL,link){
  
  resultslist$xhat1<-RL$xhat1
  resultslist$xhat2<-RL$xhat2
  resultslist$xhat3<-RL$xhat3
  resultslist$xhat4<-RL$xhat4
  
  save(resultslist,file=link)
  
}

AddSum<-function(RL){
  sumtable<-RL$fit.summary
  RL$fit.summary<-NULL
  
  RL$gompertz<-array(RL$gompertz, c(dim(RL$gompertz)[1],2,3,2))
   
  ### n_eff ###
  RL$neffM<-cbind(sumtable[grep("M_i.*,1]",row.names(sumtable)),9],
                  sumtable[grep("M_i.*,2]",row.names(sumtable)),9])
  RL$neffDelta<-cbind(sumtable[grep("Delta_i.*,1]",row.names(sumtable)),9],
                      sumtable[grep("Delta_i.*,2]",row.names(sumtable)),9])
  RL$nefftau<-cbind(sumtable[grep(paste0(variancer,".*,1]"),row.names(sumtable)),9],
                    sumtable[grep(paste0(variancer,".*,2]"),row.names(sumtable)),9],
                    sumtable[grep(paste0(variancer,".*,3]"),row.names(sumtable)),9],
                    sumtable[grep(paste0(variancer,".*,4]"),row.names(sumtable)),9])
  RL$neffbeta<-sumtable[grep("beta",row.names(sumtable)),9]
  tmp<-grep("gompertz",row.names(sumtable))
  RL$neffgompertz<-aperm(array(sumtable[tmp,9], c(2,3,2)))
  #############

  ### Rhat ###
  RL$RhatM<-cbind(sumtable[grep("M_i.*,1]",row.names(sumtable)),10],
                  sumtable[grep("M_i.*,2]",row.names(sumtable)),10])
  RL$RhatDelta<-cbind(sumtable[grep("Delta_i.*,1]",row.names(sumtable)),10],
                      sumtable[grep("Delta_i.*,2]",row.names(sumtable)),10])
  RL$Rhattau<-cbind(sumtable[grep(paste0(variancer,".*,1]"),row.names(sumtable)),10],
                    sumtable[grep(paste0(variancer,".*,2]"),row.names(sumtable)),10],
                    sumtable[grep(paste0(variancer,".*,3]"),row.names(sumtable)),10],
                    sumtable[grep(paste0(variancer,".*,4]"),row.names(sumtable)),10])
  RL$Rhatbeta<-sumtable[grep("beta",row.names(sumtable)),10]
  tmp<-grep("gompertz",row.names(sumtable))
  RL$Rhatgompertz<-aperm(array(sumtable[tmp,10], c(2,3,2)))
  #############
  
  return(RL)
}

# Calculate the centering parameters based on the posterior distribution samples from the HMC algorithm
XhatCal<-function(RL){
  
  if(RL$FRS) {
    list_nhanes<-list_nhanesFRS
  } else list_nhanes<-list_nhanesA
  
  Time<-list_nhanes$T
  age<-list_nhanes$age
  gender<-list_nhanes$male+1
  ethn<-list_nhanes$black+2*list_nhanes$white+3*list_nhanes$other
  
  FR_score<-!is.na(RL$FRSt)
  
  if (FR_score){
    
    if(RL$FRSt=="ATP") {
      FRS<-list_nhanes$FRS.ATP
    } else FRS<-list_nhanes$FRS.1998
    FRSc<-log(mean(exp(FRS)))
    
  }
  
  if(RL$eventall) {
    delta<-list_nhanes$eventall
  } else delta<-list_nhanes$eventCVDHrt
  
  # xhat<-t(cbind(RL$xhat1,RL$xhat2,RL$xhat3,RL$xhat4))
  
  beta<-RL$beta
  M_i_S<-RL$M_i_S
  M_i_D<-RL$M_i_D
  D_i_S<-RL$D_i_S
  D_i_D<-RL$D_i_D
  tau_C_S<-RL$tau_C_S
  tau_H_S<-RL$tau_H_S
  tau_C_D<-RL$tau_C_D
  tau_H_D<-RL$tau_H_D
  
  gompz<-array(RL$gompertz, c(dim(RL$gompertz)[1],2,3,2))
  
  B<-theta<-array(0,dim=dim(D_i_S))
  for (j in 1:length(gender)){
    B[,j]<-gompz[,gender[j],ethn[j],1]
    theta[,j]<-gompz[,gender[j],ethn[j],2]
  }
  rm(gompz,gender,ethn)
  
  
  # xhat is defined differently based on FRS or nhanesA(inc. FRSpop)
  if(is.null(FRS)){
    
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
  
  print("OLD XHAT = ")
  print(RL$xhat1)
  print(RL$xhat2)
  print(RL$xhat3)
  print(RL$xhat4)
  
  print("xhat = ")
  print(xhat)       
  print(" ")
  
  RL$xhat1<-xhat[1,]
  RL$xhat2<-xhat[2,]
  RL$xhat3<-xhat[3,]
  RL$xhat4<-xhat[4,]
  
  return(RL)
  
}

load(paste0(directory,"Data_cleaned/nhanes_cleaned_lists.RData"))

##########################################
################# READ-IN ################
##########################################

# sigma<-TRUE
sigma<-T
if(sigma) variancer<-"sigma" else variancer<-"tau"

# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesA_eventCVDHrt_4_3.Rdata'))
# load('Samples/savedRun_nhanesA_eventCVDHrt_0_2.Rdata')
# load('Samples/exp_betasnhanesA__RD2eventCVDHrt2.Rdata')
load(paste0(directory,'Samples/exp_betasnhanesA_eventCVDHrt_sigma_4.Rdata'))
R1<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL1<-resultslist
RL1$FRS<-FALSE
RL1$eventall<-FALSE
RL1$pop<-FALSE
RL1$FRSt<-NA
RL1<-AddSum(RL1)
RL1$fit.summary<-NULL
# RL1%<>%XhatCal()
# saveresultslist(resultslist,RL2,paste0(directory,'Samples/exp_betasnhanesA_eventCVDHrt10.Rdata'))
# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesA_eventall_3_1.Rdata'))
load(paste0(directory,'Samples/exp_betasnhanesA_eventall_sigma_4.Rdata'))
R2<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL2<-resultslist
RL2$FRS<-FALSE
RL2$eventall<-TRUE
RL2$pop<-FALSE
RL2$FRSt<-NA
RL2<-AddSum(RL2)
RL2$fit.summary<-NULL
# RL2%<>%XhatCal()
# saveresultslist(resultslist,RL1,paste0(directory,'Samples/exp_betasnhanesA_eventall10.Rdata'))
# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesFRSpop_eventCVDHrt_4_3.Rdata'))
load(paste0(directory,'Samples/exp_betasnhanesFRSpop_eventCVDHrt_sigma_4.Rdata'))
R3<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL3<-resultslist
RL3$FRS<-TRUE
RL3$eventall<-FALSE
RL3$pop<-TRUE
RL3$FRSt<-NA
RL3<-AddSum(RL3)
RL3$fit.summary<-NULL
# RL3%<>%XhatCal()
# saveresultslist(resultslist,RL4,paste0(directory,'Samples/exp_betasnhanesFRSpop_eventCVDHrt10.Rdata'))
# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesFRSpop_eventall_4_3.Rdata'))
load(paste0(directory,'Samples/exp_betasnhanesFRSpop_eventall_sigma_4.Rdata'))
R4<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL4<-resultslist
RL4$FRS<-TRUE
RL4$eventall<-TRUE
RL4$pop<-TRUE
RL4$FRSt<-NA
RL4<-AddSum(RL4)
RL4$fit.summary<-NULL
# RL4%<>%XhatCal()
# saveresultslist(resultslist,RL3,paste0(directory,'Samples/exp_betasnhanesFRSpop_eventall10.Rdata'))
# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesFRS_eventCVDHrt_FRSATP_4_2.Rdata'))
# load('/media/patten/My Passport/BP_Nhanes/Samples/savedRun_nhanesFRS_eventCVDHrt_FRSATP_0_1.Rdata')
# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesFRS_eventCVDHrt_FRSATP_0_1.Rdata'))
load(paste0(directory,'Samples/exp_betasnhanesFRS_eventCVDHrt_FRSATP_sigma_6.Rdata'))
R5<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL5<-resultslist
RL5$FRS<-TRUE
RL5$eventall<-FALSE
RL5$pop<-FALSE
RL5$FRSt<-'ATP'
RL5<-AddSum(RL5)
RL5$fit.summary<-NULL
# RL5%<>%XhatCal()
# saveresultslist(resultslist,RL8,paste0(directory,'Samples/exp_betasnhanesFRS_eventCVDHrt_FRSATP10.Rdata'))
# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesFRS_eventall_FRSATP_3_1.Rdata'))
load(paste0(directory,'Samples/exp_betasnhanesFRS_eventall_FRSATP_sigma_4.Rdata'))
R6<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL6<-resultslist
RL6$FRS<-TRUE
RL6$eventall<-TRUE
RL6$pop<-FALSE
RL6$FRSt<-'ATP'
RL6<-AddSum(RL6)
RL6$fit.summary<-NULL
# RL6%<>%XhatCal()
# saveresultslist(resultslist,RL7,paste0(directory,'Samples/exp_betasnhanesFRS_eventall_FRSATP10.Rdata'))
# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesFRS_eventCVDHrt_FRS1998_4_4.Rdata'))
load(paste0(directory,'Samples/exp_betasnhanesFRS_eventCVDHrt_FRS1998_sigma_6.Rdata'))
R7<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL7<-resultslist
RL7$FRS<-TRUE
RL7$eventall<-FALSE
RL7$pop<-FALSE
RL7$FRSt<-'1998'
RL7<-AddSum(RL7)
RL7$fit.summary<-NULL
# RL7%<>%XhatCal()
# saveresultslist(resultslist,RL6,paste0(directory,'Samples/exp_betasnhanesFRS_eventCVDHrt_FRS199810.Rdata'))
# load(paste0(directory,'Samples/summariesdata/savedRun_nhanesFRS_eventall_FRS1998_3_4.Rdata'))
load(paste0(directory,'Samples/exp_betasnhanesFRS_eventall_FRS1998_sigma_4.Rdata'))
R8<-t(cbind(resultslist$xhat1,resultslist$xhat2,resultslist$xhat3,resultslist$xhat4))
RL8<-resultslist
RL8$FRS<-TRUE
RL8$eventall<-TRUE
RL8$pop<-FALSE
RL8$FRSt<-'1998'
RL8<-AddSum(RL8)
RL8$fit.summary<-NULL
# RL8%<>%XhatCal()
# saveresultslist(resultslist,RL5,paste0(directory,'Samples/exp_betasnhanesFRS_eventall_FRS199810.Rdata'))
rm(resultslist,AddSum)

if(sigma){
  tity<-c(a="M Home-Clinic Mean",b="$\\Delta$ Home-Clinic Difference",c="$\\tau$ Clinic Precision",d="$\\tau$ Home Precision")
  texty<-c("M","Delta","tauC","tauH")
  titles<-c("beta[1] - M Systolic or FRS value","beta[2] - Delta Systolic","beta[3] - tau Clinic Systolic","beta[4] - tau Home Systolic",
            "beta[5] - M Diastolic","beta[6] - Delta Diastolic","beta[7] - tau Clinic Diastolic","beta[8] - tau Home Diastolic")
  titlesL<-c("$\\beta_1$ - M Systolic or FRS value","$\\beta_2$ - $\\Delta$ Systolic","$\\beta_3$ - $\\tau$ Clinic Systolic","$\\beta_4$ - $\\tau$ Home Systolic",
             "$\\beta_5$ - M Diastolic","$\\beta_6$ - $\\Delta$ Diastolic","$\\beta_7$ - $\\tau$ Clinic Diastolic","$\\beta_8$ - $\\tau$ Home Diastolic")
} else {
  tity<-c(a="M Home-Clinic Mean",b="$\\Delta$ Home-Clinic Difference",c="$\\sigma$ Clinic Std. Dev.",d="$\\sigma$ Home Std. Dev.")
  texty<-c("M","Delta","sigmaC","sigmaH")
  titles<-c("beta[1] - M Systolic or FRS value","beta[2] - Delta Systolic","beta[3] - sigma Clinic Systolic","beta[4] - sigma Home Systolic",
            "beta[5] - M Diastolic","beta[6] - Delta Diastolic","beta[7] - sigma Clinic Diastolic","beta[8] - sigma Home Diastolic")
  titlesL<-c("$\\beta_1$ - M Systolic or FRS value","$\\beta_2$ - $\\Delta$ Systolic","$\\beta_3$ - $\\sigma$ Clinic Systolic","$\\beta_4$ - $\\sigma$ Home Systolic",
             "$\\beta_5$ - M Diastolic","$\\beta_6$ - $\\Delta$ Diastolic","$\\beta_7$ - $\\sigma$ Clinic Diastolic","$\\beta_8$ - $\\sigma$ Home Diastolic")
}

##########################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%% COMPARE DIFFERENT SIMULATIONS (DIFFERENT INPUT DATA) USING PLOTS %%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
############### BETA PLOTS ###############

beta<-data.frame()
for(i in 1:8){
  RL<-get(paste0("RL",i))
  if(is.na(RL$FRSt)) btmp<-ConvertBeta(RL) else {
    if(RL$FRSt=="ATP") FRS<-list_nhanesFRS$FRS.ATP else FRS<-list_nhanesFRS$FRS.1998
    btmp<-ConvertBeta(RL,FRS=FRS)
  }
  beta%<>%rbind(cbind(btmp,data.frame(runnum=rep(i,nrow(btmp)))))
}
rm(btmp,FRS,RL)
beta$runnum%<>%as.factor()

beta_sum<-beta%>%group_by(runnum,variable)%>%summarise(nmean=mean(normalised),
                                             nsd=sd(normalised),
                                             nhp95=quantile(normalised,0.95),
                                             nhp05=quantile(normalised,0.05),
                                             vmean=mean(value),
                                             vsd=sd(value),
                                             vhp95=quantile(value,0.95),
                                             vhp05=quantile(value,0.05))

names(beta_sum)<-c("Run Number","Variable","Mean (Normalised)","SD (Normalised)","HP95 (Normalised)","HP05 (Normalised)","Mean","SD","HP95","HP05")
write_csv(beta_sum,file="Results/beta_all.csv")

p<-ggplot(data=beta,aes(x=runnum,y=normalised,fill=variable))+geom_violin()+geom_abline(slope=0,intercept = 0) +#+ggtitle()
  theme(plot.title = element_text(hjust = 0.5))+ylab(TeX("$\\beta$ Linear Predictor Parameter")) +
  xlab("Run Number") + ggtitle(TeX("Normalised $\\beta$ Parameter"))
p<-p+facet_wrap( ~ variable, nrow = 4, scales = "free")+ theme(plot.title = element_text(hjust = 0.5));p
ggsave('Beta_parameter_normalised.png', plot=p,path = paste0(directory,'Plots/beta'),width = 15,height = 9)

gompertz<-data.frame()
for(i in 1:8){
  RL<-get(paste0("RL",i))
  gompz<-ConvertGompertz(RL$gompertz)
  gompertz%<>%rbind(cbind(gompz,data.frame(runnum=rep(i,nrow(gompz)))))
}
rm(gompz)
gompertz$runnum%<>%as.factor()

p<-ggplot(data=gompertz,aes(x=Ethnicity,y=B,fill=Gender))+geom_violin()+scale_y_log10()+#+ggtitle()
  theme(plot.title = element_text(hjust = 0.5))+ylab("B - Gompertz Parameter") + xlab("Ethnicity")
p<-p+facet_wrap( ~ runnum, nrow = 4, scales = "fixed")+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Run Number");p
ggsave('B_parameter.png', plot=p,path = paste0(directory,'Plots/gompertz'),width = 13,height = 9)

p<-ggplot(data=gompertz,aes(x=Ethnicity,y=theta,fill=Gender))+geom_violin()+#+ggtitle()
  theme(plot.title = element_text(hjust = 0.5))+ylab(TeX("$\\theta$ - Gompertz Parameter")) + xlab("Ethnicity")
p<-p+facet_wrap( ~ runnum, nrow = 4, scales = "fixed")+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Run Number");p
ggsave('theta_parameter.png', plot=p,path = paste0(directory,'Plots/gompertz'),width = 13,height = 9)

############### XHAT PLOTS ###############
xhat<-vector()
DS<-vector()
runnum<-vector()
indy<-vector()
event<-vector()
xhatD<-xhatS<-numeric(8)

sss<-c(a="$\\hat{x}_1$ - M or FRS value",b="$\\hat{x}_2$ - $\\Delta$",c="$\\hat{x}_3$ - $\\tau$ Clinic",d="$\\hat{x}_4$ - $\\tau$ Home")
for (j in 1:4){
  strng=sss[j]
  for(i in 1:8){
    RL<-get(paste0("RL",i))
    if(j==1 & !is.na(RL$FRSt)){
      xhatS[i]<-NA
      xhatD[i]<-NA
    } else{
      xhatS[i]<-RL[[paste0("xhat",j)]][1]
      xhatD[i]<-RL[[paste0("xhat",j)]][2]
    }    
  }  
  xhat<-c(xhat,xhatS,xhatD)
  DS<-c(DS,rep("Systolic",8),rep("Diastolic",8))
  runnum<-c(runnum,1:8,1:8)
  indy<-c(indy,rep(strng,16))
  event<-c(event,rep(c("All","CVDHrt"),8))
}
xhatty<-data.frame(xhat=xhat,DS=DS,runnum=as.factor(runnum),indy=indy,event=event)
rm(RL,xhatS,xhatD,runnum,DS,xhat,event,indy)

p<-ggplot(xhatty,aes(x=runnum,y=xhat,fill=DS,color=DS))+geom_point(aes(shape=event,stroke=3))+ggtitle(TeX(titlesL[j]))+ylab(TeX("$\\hat{x}$ Mean Adjustment"))+xlab("Run Number");
p<-p+facet_wrap( ~ DS+indy, scales = "free_y",labeller = as_labeller(TeX,default = label_parsed), nrow = 2)+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle(TeX("$\\hat{x}$ Index"));p
ggsave('xhat_all.png', plot=p,path = paste0(directory,'Plots/xhat'),width = 13,height = 7)

# png(paste0(directory,'Plots/xhat/M_hat.png'))
# plot(c(RL1$xhat1,RL2$xhat1,RL3$xhat1,RL4$xhat1,RL5$xhat1,RL6$xhat1,RL7$xhat1,RL8$xhat1)[seq(1,16,2)],pch=21,type='p',col='red',ylab = 'M-hat Systolic',xlab="Run Number")
# points(c(RL1$xhat1,RL2$xhat1,RL3$xhat1,RL4$xhat1,RL5$xhat1,RL6$xhat1,RL7$xhat1,RL8$xhat1)[seq(2,16,2)],pch=22,col='blue',ylab = 'Mhat Diastolic')
# legend("topright", legend = c("Systolic", "Diastolic"),col=c('red','blue'),pch=c(21,22))
# dev.off()

############### N_EFF VIOLIN PLOTS ###############
variables<-c("neffM","neffDelta")
j=1
for (xxx in variables){
  LLL<-length(RL1[[xxx]][,1])
  DF<-data.frame(M=c(RL1[[xxx]][,1],RL1[[xxx]][,2]),DS=c(rep("Systolic",LLL),rep("Diastolic",LLL)),runnum=rep(1,2*LLL),FRS=rep(RL1$FRS,2*LLL))
  for(i in 2:8){
    RL<-get(paste0("RL",i))
    LLL<-length(RL[[xxx]][,1])
    DF<-rbind(DF,data.frame(M=c(RL[[xxx]][,1],RL[[xxx]][,2]),DS=c(rep("Systolic",LLL),rep("Diastolic",LLL)),runnum=rep(i,2*LLL),FRS=rep(RL$FRS,2*LLL)))
  }
  p<-ggplot(data=DF,aes(x=as.factor(runnum),y=M,fill=DS))+geom_violin()+xlab("Run Number")+ggtitle(TeX(tity[j]))+theme(plot.title = element_text(hjust = 0.5))+ylab("N_eff")
  ggsave(paste0('N-eff_',texty[j],'.png'), plot=p,path = paste0(directory,'Plots/neff'),width = 10,height = 3.)
  j=j+1
}

### N_EFF BETA ###
#x- runnum, y-n_eff, multi plots for each beta value - 4 columns, 2 rows, split by DS and param
variables<-c("neffbeta","Rhatbeta")
DS<-c(rep("Systolic",4),rep("Diastolic",4))
param<-c(a="M",b="$\\Delta$",c=paste0("$\\",variancer,"$ Clinic"),d=paste0("$\\",variancer,"$ Home")); param<-c(param,param)
DF<-data.frame(M=RL1[[variables[1]]],N=RL1[[variables[2]]],DS=DS,param=param,runnum=rep(1,8),FRS=rep(RL1$FRS,8))
for(i in 2:8){
  RL<-get(paste0("RL",i))
  if (!is.na(RL$FRSt)){tmp1<-c(RL[[variables[1]]][1:4],NA,RL[[variables[1]]][5:7]);tmp2<-c(RL[[variables[2]]][1:4],NA,RL[[variables[2]]][5:7])}
  else {tmp1<-RL[[variables[1]]];tmp2<-RL[[variables[2]]]}
  DF<-rbind(DF,data.frame(M=tmp1,N=tmp2,DS=DS,param=param,runnum=rep(i,8),FRS=rep(RL$FRS,8)))
}
p<-ggplot(data=DF,aes(x=as.factor(runnum),y=M))+geom_point(aes(colour=param,shape=DS))+xlab("Run Number")+ggtitle(TeX(param))+theme(plot.title = element_text(hjust = 0.5))+ylab("N_eff")
p<-p+facet_wrap( ~ param+DS, labeller = as_labeller(TeX,default = label_parsed), nrow = 1)+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("N_eff");p
ggsave('N-eff_beta.png', plot=p,path = paste0(directory,'Plots/neff'),width = 13,height = 3.)
p<-ggplot(data=DF,aes(x=as.factor(runnum),y=N,fill=DS))+geom_point(aes(colour=param,shape=DS))+xlab("Run Number")+ggtitle(TeX(param))+theme(plot.title = element_text(hjust = 0.5))+ylab("R-hat")
p<-p+facet_wrap( ~ param+DS, labeller = as_labeller(TeX,default = label_parsed), nrow = 1)+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("R-hat");p
ggsave('Rhat_beta.png', plot=p,path = paste0(directory,'Plots/neff'),width = 13,height = 3.)

### N_EFF GOMPERTZ ###
#x- runnum, y-n_eff, multi plots for each beta value - 4 columns, 2 rows, split by DS and param
variables<-c("neffgompertz","Rhatgompertz")
gender<-rep(c(rep("Female",3),rep("Male",3)),2)
ethn<-rep(c("Black","White","Other"),4)
DS<-c(rep("B",6),rep("$\\theta$",6))
DF<-data.frame(M=unlist(list(RL1[[variables[1]]])),N=unlist(list(RL1[[variables[2]]])),parameter=DS,gender=gender,ethn=ethn,runnum=rep(1,12))
for(i in 2:8){
  RL<-get(paste0("RL",i))
  DF<-rbind(DF,data.frame(M=unlist(list(RL[[variables[1]]])),N=unlist(list(RL[[variables[2]]])),parameter=DS,gender=gender,ethn=ethn,runnum=rep(i,12)))
}
p<-ggplot(data=DF,aes(x=as.factor(runnum),y=M,fill=parameter))+geom_point(aes(colour=parameter))+xlab("Run Number")+ggtitle(TeX(param))+theme(plot.title = element_text(hjust = 0.5))+ylab("N_eff")
p<-p+facet_wrap( ~ gender+ethn+parameter, labeller = as_labeller(TeX,default = label_parsed), nrow = 2)+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("N_eff");p
ggsave('N-eff_gompertz.png', plot=p,path = paste0(directory,'Plots/neff'),width = 13,height = 7.)
p<-ggplot(data=DF,aes(x=as.factor(runnum),y=N,fill=parameter))+geom_point(aes(colour=parameter))+xlab("Run Number")+ggtitle(TeX(param))+theme(plot.title = element_text(hjust = 0.5))+ylab("R-hat")
p<-p+facet_wrap( ~ gender+ethn+parameter, labeller = as_labeller(TeX,default = label_parsed), nrow = 2)+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("R-hat");p
ggsave('Rhat_gompertz.png', plot=p,path = paste0(directory,'Plots/neff'),width = 13,height = 7.)

#######################   IMPORT DATA   ##############################
# load(paste0(directory,"Data_cleaned/nhanes_cleaned_frames.RData"))
xSC_NF<-list_nhanesA$sys
xDC_NF<-list_nhanesA$dias
xSH_NF<-list_nhanesA$sys_home
xDH_NF<-list_nhanesA$dias_home
xSC_F<-list_nhanesFRS$sys
xDC_F<-list_nhanesFRS$dias
xSH_F<-list_nhanesFRS$sys_home
xDH_F<-list_nhanesFRS$dias_home
### Recalculate initial xhat values ###
M_av_NF <- c(mean(rowMeans( as.matrix(xSC_NF) )),mean(rowMeans( as.matrix(xDC_NF) )))
Delta_av_NF <- c(mean( rowMeans( abs((as.matrix(xSC_NF)-as.matrix(xSH_NF))/2 ) )),mean( rowMeans( abs((as.matrix(xDC_NF)-as.matrix(xDH_NF))/2 ) )))
sds1 <- max(apply(xSC_NF, 1, sd))
sds2 <- max(apply(xSH_NF, 1, sd))
sds3 <- max(apply(xDC_NF, 1, sd))
sds4 <- max(apply(xDH_NF, 1, sd))
tauis_av_NF <- c(mean(1/sds1^2) , mean(1/sds2^2), mean(1/sds3^2), mean(1/sds4^2))
tauis_NF <- c(1/sds1^2 , 1/sds2^2, 1/sds3^2, 1/sds4^2)
# FRS
M_av_F <- c(mean(rowMeans( as.matrix(xSC_F) )),mean(rowMeans( as.matrix(xDC_F) )))
Delta_av_F <- c(mean( rowMeans( abs((as.matrix(xSC_F)-as.matrix(xSH_F))/2 ) )),mean( rowMeans( abs((as.matrix(xDC_F)-as.matrix(xDH_F))/2 ) )))
sds1 <- max(apply(xSC_F, 1, sd))
sds2 <- max(apply(xSH_F, 1, sd))
sds3 <- max(apply(xDC_F, 1, sd))
sds4 <- max(apply(xDH_F, 1, sd))
tauis_av_F <- c(mean(1/sds1^2) , mean(1/sds2^2), mean(1/sds3^2), mean(1/sds4^2))
tauis_F <- c(1/sds1^2 , 1/sds2^2, 1/sds3^2, 1/sds4^2)
rm(sds1,sds2,sds3,sds4)
###################################################################

Nits<-5

xhatS<-xhatD<-array(0,dim = c(Nits,4))
DF<-data.frame()
### Plot Convergence of xhat values ###
for(i in 1:8){
  
  namer<-paste0("RL",i)
  RL<-get(namer)
  
  if (RL$FRS){FF<-"FRS"} else {FF<-"A"}
  if (RL$pop){popy<-"pop"} else {popy<-""}
  if (RL$eventall){ev<-"all"} else {ev<-"CVDHrt"}
  if (!is.na(RL$FRSt)){add<-paste0("_FRS",RL$FRSt)} else {add<-""}
  
  # Run Number 1 (eventall & Nhanes):
  if(!is.na(RL$FRSt)){
    xhatS[1,1]<-M_av_F[1];      xhatD[1,1]<-M_av_F[2]
    xhatS[1,2]<-Delta_av_F[1];  xhatD[1,2]<-Delta_av_F[2]
    xhatS[1,3]<-tauis_av_F[1];  xhatD[1,3]<-tauis_av_F[3]
    xhatS[1,4]<-tauis_av_F[2];  xhatD[1,4]<-tauis_av_F[4] 
  } else {
    xhatS[1,1]<-M_av_NF[1];      xhatD[1,1]<-M_av_NF[2]
    xhatS[1,2]<-Delta_av_NF[1];  xhatD[1,2]<-Delta_av_NF[2]
    xhatS[1,3]<-tauis_av_NF[1];  xhatD[1,3]<-tauis_av_NF[3]
    xhatS[1,4]<-tauis_av_NF[2];  xhatD[1,4]<-tauis_av_NF[4]
  }
  
  for (j in 1:(Nits-1)){
    print(paste0(directory,'Samples/exp_betasnhanes',FF,popy,'_event',ev,add,j,'.Rdata'))
    tryCatch(load(paste0(directory,'Samples/exp_betasnhanes',FF,popy,'_event',ev,add,j,'.Rdata')),error = function(e) resultslist<-NULL)
    if(is.null(resultslist)) {
      xhatS[j+1,]<-xhatD[j+1,]<-NA
    } else {
      resultslist$xhat1<-as.numeric(resultslist$xhat1)
      resultslist$xhat2<-as.numeric(resultslist$xhat2)
      resultslist$xhat3<-as.numeric(resultslist$xhat3)
      resultslist$xhat4<-as.numeric(resultslist$xhat4)
      xhatS[j+1,]<-c(resultslist$xhat1[1],resultslist$xhat2[1],resultslist$xhat3[1],resultslist$xhat4[1])
      xhatD[j+1,]<-c(resultslist$xhat1[2],resultslist$xhat2[2],resultslist$xhat3[2],resultslist$xhat4[2])
    }
  }
  if (RL$FRS&!RL$pop){xhatS[,1]<-xhatD[,1]<-NA}
  
  DF<-rbind(DF,data.frame(xhat=c(unlist(list(t(xhatS))),unlist(list(t(xhatD)))),DS=c(rep("Systolic",Nits*4),rep("Diastolic",Nits*4)),parameter=rep(tity,Nits),
                          iteration=c(rep(1,Nits),rep(2,Nits),rep(3,Nits),rep(4,Nits)),simulation=as.factor(rep(namer,Nits*8)),event=as.factor(rep(ev,Nits*8))))
}

p<-ggplot(data=DF,aes(x=as.factor(iteration), y=xhat, colour=simulation, shape=event))+geom_point()+xlab("Iteration Number")+theme(plot.title = element_text(hjust = 0.5))+ylab(TeX("$\\hat{x}$ Value"))
p<-p+facet_wrap( ~ DS+parameter, labeller = as_labeller(TeX,default = label_parsed), nrow = 2, scales = "free_y")+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle(TeX("Convergence in $\\hat{x}$ Nhanes Value"));p
ggsave('xhat_conv.png', plot=p,path = paste0(directory,'Plots/xhat'),width = 10,height = 5.)

rm(xhatS,xhatD,DF)

############################### ROC #################################'

survy<-GetSurvival(RL = RL,roc = F,Gender="female",usexhat=TRUE)
 
GetPredictions<-function(roc=F,Ethnicity=NULL,Gender=NULL){
  ### DATA ROC ###
  preds<-data.frame()
  for(i in 1:8){
    RL<-get(paste0("RL",i))
    survy<-GetSurvival(RL = RL,roc = roc,Ethnicity=Ethnicity,Gender=Gender,usexhat=TRUE)
    preds<-rbind(preds,cbind(survy,runnum=rep(i,dim(survy)[1])))
  }
  
  preds$runnum<-as.factor(preds$runnum)
  return(preds)
  
}

# preds<-GetPredictions(F)

GetPredsEthnicity<-function(roc=F){
  
  EEE<-preds<-data.frame()
  for (Ethnicity in c("black","white","other")){
    preds%<>%rbind(GetPredictions(roc=roc,Ethnicity=Ethnicity))
    EEE<-rbind(EEE,data.frame(ethn=rep(Ethnicity,(nrow(preds)-nrow(EEE)))))
  }
  
  preds%>%cbind(EEE)
  
}
 
Epreds<-GetPredsEthnicity()

GetPredsGender<-function(roc=F){
  
  GGG<-preds<-data.frame()
  for (Gender in c("female","male")){
    preds%<>%rbind(GetPredictions(roc=roc,Gender=Gender))
    GGG<-rbind(GGG,data.frame(genre=rep(Gender,(nrow(preds)-nrow(GGG)))))
  }
  
  preds%>%cbind(GGG)
  
}
 
Gpreds<-GetPredsGender()

preds<-GetPredictions(TRUE)

preds$survived<-preds$d_j
preds$pred<-preds$H_j

for (i in as.integer(unique(preds$runnum))){
  tmp<-dplyr::select(filter(preds,runnum==i),-runnum)
  
  roc<-calculate_roc(df=tmp,cost_of_fp =1,cost_of_fn=1,n = 300)
  p<-plot_roc(roc,0.75,1,1)
  ggsave(paste0("redlinpred_ROC-RL",i,'.png'), plot=p,path = paste0(directory,'Plots/Survival'),width = 10,height = 5.)
  
  p1<-ggplot(tmp,aes(log(pred),group=survived,fill=as.factor(survived)))+geom_density(alpha=0.3) + #scale_x_continuous(trans="log10") +
    scale_fill_discrete(name = "Mortality", labels = c("Alive","Dead")) + xlab("log(H) - Hazard") + ylab("Density")
  # p2<-ggplot(tmp,aes(pred))+geom_density(alpha=0.3) + scale_x_continuous(trans="log10") +
  #   scale_fill_discrete(name = "Mortality", labels = c("Alive","Dead")) + xlab("H - Hazard") + ylab("Density")
  # p<-ggarrange(p1,p2,nrow=1)
  ggsave(paste0("redlinpred_Hist_survival-RL",i,".png"), plot=p1,path = paste0(directory,'Plots/Survival'),width = 7,height = 5.)
  # p1
  
}

preds<-GetPredictions(F)

# preds%>%group_by(runnum)%>%summarise(RMSE=sqrt(mean((H_j-d_j)*(H_j-d_j)/d_j,na.rm = T)))

p<-ggplot(preds,aes(x=H_j,y=d_j))+geom_line() + geom_point() + #scale_fill_discrete(name = "Gender", labels = c("Male","Female"))+
  xlab("Cumulative Hazard") + ylab("Cumulative Deaths") + #scale_color_manual(values = c(1="Male",2="Female"))
  geom_abline(intercept=0, slope=1,colour="black");
p<-p+facet_wrap( ~ runnum, nrow = 2, scales = "free")+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Run Number");p
ggsave("redlinpred_Cumulative_haz-death_age.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 10,height = 5.)

Epreds$Ethnicity<-Epreds$ethn ; Epreds$ethn<-NULL
p<-ggplot(Epreds,aes(x=H_j,y=d_j,group=Ethnicity))+geom_point(aes(color=Ethnicity),size=1) + geom_line(aes(color=Ethnicity)) +#scale_fill_discrete(name = "Gender", labels = c("Male","Female"))+
  xlab("Cumulative Hazards") + ylab("Cumulative Events") + #scale_color_manual(values = c(1="Male",2="Female"))
  geom_abline(intercept=0, slope=1,colour="black")
p<-p+facet_wrap( ~ runnum, nrow = 2, scales = "free")+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Run Number");p
ggsave("Cumulative_ev-haz_age-srt_ethnicity.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 10,height = 5.)

Gpreds$Gender<-Gpreds$genre ; Gpreds$genre<-NULL
p<-ggplot(Gpreds,aes(x=H_j,y=d_j,group=Gender))+geom_point(aes(color=Gender),size=1) + geom_line(aes(color=Gender)) +#scale_fill_discrete(name = "Gender", labels = c("Male","Female"))+
  xlab("Cumulative Hazards") + ylab("Cumulative Events") + #scale_color_manual(values = c(1="Male",2="Female"))
  geom_abline(intercept=0, slope=1,colour="black")
p<-p+facet_wrap( ~ runnum, nrow = 2, scales = "free")+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Run Number");p
ggsave("Cumulative_ev-haz_age-srt_gender.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 10,height = 5.)


#################### BAYES FACTORS #######################
BF<-data.frame()
for(i in 1:8){
  RL<-get(paste0("RL",i))
  BF%<>%rbind(cbind(RL$BF,data.frame(runnum=rep(i,length(RL$BF$H_1)))))
}

p<-ggplot(BF,aes(x=X1,y=BF))+geom_point(aes(colour=H_1))+xlab(TeX("$\\beta$ Index"))+ylab("Bayes Factor")#+scale_fill_discrete(name = "BF", labels = c("p<0","p>0","No Value"));
p<-p+facet_wrap( ~ runnum, nrow = 2, scales = "free_y")+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Run Number");p

p<-ggplot(BF,aes(x=runnum,y=BF))+geom_point(aes(colour=H_1))+xlab("Run Number")+ylab("Bayes Factor")#+scale_fill_discrete(name = "BF", labels = c("p<0","p>0","No Value"));p
p<-p+facet_wrap( ~ BF$, nrow = 2, scales = "fixed")+scale_y_log10()+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle(TeX("$\\beta$ Index"));p
ggsave("BF_beta.png", plot=p,path = paste0(directory,'Plots/BF'),width = 12,height = 5.)


#%%% Lexis diagrams of the population %%%#
plot_discrete_cbar = function(
  breaks, # Vector of breaks. If +-Inf are used, triangles will be added to the sides of the color bar
  palette = "Greys", # RColorBrewer palette to use
  colors = RColorBrewer::brewer.pal(length(breaks) - 1, palette), # Alternatively, manually set colors
  direction = 1, # Flip colors? Can be 1 or -1
  spacing = "natural", # Spacing between labels. Can be "natural" or "constant"
  border_color = NA, # NA = no border color
  legend_title = NULL,
  legend_direction = "horizontal", # Can be "horizontal" or "vertical"
  font_size = 5,
  expand_size = 1, # Controls spacing around legend plot
  spacing_scaling = 1, # Multiplicative factor for label and legend title spacing
  width = 0.2, # Thickness of color bar
  triangle_size = 0.000001 # Relative width of +-Inf triangles
) {
  require(ggplot2)
  if (!(spacing %in% c("natural", "constant"))) stop("spacing must be either 'natural' or 'constant'")
  if (!(direction %in% c(1, -1))) stop("direction must be either 1 or -1")
  if (!(legend_direction %in% c("horizontal", "vertical"))) stop("legend_direction must be either 'horizontal' or 'vertical'")
  breaks = as.numeric(breaks)
  new_breaks = sort(unique(breaks))
  if (any(new_breaks != breaks)) warning("Wrong order or duplicated breaks")
  breaks = new_breaks
  if (class(colors) == "function") colors = colors(length(breaks) - 1)
  if (length(colors) != length(breaks) - 1) stop("Number of colors (", length(colors), ") must be equal to number of breaks (", length(breaks), ") minus 1")
  if (!missing(colors)) warning("Ignoring RColorBrewer palette '", palette, "', since colors were passed manually")
  
  if (direction == -1) colors = rev(colors)
  
  inf_breaks = which(is.infinite(breaks))
  if (length(inf_breaks) != 0) breaks = breaks[-inf_breaks]
  plotcolors = colors
  
  n_breaks = length(breaks)
  
  labels = breaks
  
  if (spacing == "constant") {
    breaks = 1:n_breaks
  }
  
  r_breaks = range(breaks)
  
  cbar_df = data.frame(stringsAsFactors = FALSE,
                       y = breaks,
                       yend = c(breaks[-1], NA),
                       color = as.character(1:n_breaks)
  )[-n_breaks,]
  
  xmin = 1 - width/2
  xmax = 1 + width/2
  
  cbar_plot = ggplot(cbar_df, aes(xmin=xmin, xmax = xmax, ymin = y, ymax = yend, fill = factor(color, levels = 1:length(colors)))) +
    geom_rect(show.legend = FALSE,
              color=border_color)
  
  if (any(inf_breaks == 1)) { # Add < arrow for -Inf
    firstv = breaks[1]
    polystart = data.frame(
      x = c(xmin, xmax, 1),
      y = c(rep(firstv, 2), firstv - diff(r_breaks) * triangle_size)
    )
    plotcolors = plotcolors[-1]
    cbar_plot = cbar_plot +
      geom_polygon(data=polystart, aes(x=x, y=y),
                   show.legend = FALSE,
                   inherit.aes = FALSE,
                   fill = colors[1],
                   color=border_color)
  }
  if (any(inf_breaks > 1)) { # Add > arrow for +Inf
    lastv = breaks[n_breaks]
    polyend = data.frame(
      x = c(xmin, xmax, 1),
      y = c(rep(lastv, 2), lastv + diff(r_breaks) * triangle_size)
    )
    plotcolors = plotcolors[-length(plotcolors)]
    cbar_plot = cbar_plot +
      geom_polygon(data=polyend, aes(x=x, y=y),
                   show.legend = FALSE,
                   inherit.aes = FALSE,
                   fill = colors[length(colors)],
                   color=border_color)
  }
  
  if (legend_direction == "horizontal") { #horizontal legend
    mul = 1
    x = xmin
    xend = xmax
    cbar_plot = cbar_plot + coord_flip()
    angle = 0
    legend_position = xmax + 0.1 * spacing_scaling
  } else { # vertical legend
    mul = -1
    x = xmax
    xend = xmin
    angle = -90
    legend_position = xmax + 0.2 * spacing_scaling
  }
  
  cbar_plot = cbar_plot +
    geom_segment(data=data.frame(y = breaks, yend = breaks),
                 aes(y=y+y_space, yend=yend),
                 x = x - 0.05 * mul * spacing_scaling, xend = xend,
                 inherit.aes = FALSE) +
    annotate(geom = 'text', x = x - 0.1 * mul * spacing_scaling, y = breaks,
             label = labels,
             size = font_size) +
    scale_x_continuous(expand = c(expand_size,expand_size)) +
    scale_fill_manual(values=plotcolors) +
    theme_void()
  
  if (!is.null(legend_title)) { # Add legend title
    cbar_plot = cbar_plot +
      annotate(geom = 'text', x = legend_position, y = mean(r_breaks),
               label = legend_title,
               angle = angle,
               size = font_size)
  }
  
  cbar_plot
}
# Function to prepare the data to be in the correct format to plot later
PrepLexisData<-function(input,funcy=NULL){
  output<-data.frame()
  for(i in 1:nrow(input)){
    
    tmp<-input[i,]
    # Bottom-left
    output%<>%rbind(cbind(tmp,data.frame(id=i)))
    # Bottom-right
    tmp$sheep.yr<-tmp$sheep.yr+1
    output%<>%rbind(cbind(tmp,data.frame(id=i)))
    # Top
    tmp$age<-tmp$age+1
    output%<>%rbind(cbind(tmp,data.frame(id=i)))
    
    # Then the second triangle component
    # Top-left
    output%<>%rbind(cbind(tmp,data.frame(id=-i)))
    # Top-right
    tmp$sheep.yr<-tmp$sheep.yr+1
    output%<>%rbind(cbind(tmp,data.frame(id=-i)))
    # Bottom-left
    tmp$age<-tmp$age-1
    tmp$sheep.yr<-tmp$sheep.yr-1
    output%<>%rbind(cbind(tmp,data.frame(id=-i)))
    
  }
  # Convert from year to an actual date (it doesn't matter if it is 1st Jan or any other day)
  output$date<-paste0(as.character(output$sheep.yr),"-01-01")
  
  # Create a colour scheme for the number of births and deaths
  n<-20
  if(!is.null(funcy)) toty<-funcy(output$total) else toty<-output$total
  tmp <- hist(toty,breaks = n,plot = F) ; tmp<-findInterval(toty, tmp$breaks)
  cols<-RColorBrewer::brewer.pal(n = 11, name = "Spectral") ; cols<-colorRampPalette(cols)(n)
  output$cols<-cols[tmp]
  
  return(output)
}

mylexis <- lexis_grid(year_start = ceiling(min(list_nhanesA$T)), 
                      year_end =  ceiling(max(list_nhanesA$T)+2), 
                      age_start = 50,
                      age_end = 99)

survTot<-GetSurvival(RL=RL1,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")

p<-plot_roc(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),0.75,costy = F) +
  theme(plot.title = element_text(hjust = 0.5))+ggtitle(TeX("All $\\beta$ Linear Predictor Terms"))
ggsave("ROCAll.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 5,height = 4) 

survFrame<-data.frame()

survTot<-GetSurvival(RL=RL1,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")
tmp<-RL1; tmp$beta[]<-0
survDemog<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDemog)[1:2]<-c("pred","survived")
tmp<-RL1; tmp$beta[,2:8]<-0
survFRS<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRS)[1:2]<-c("pred","survived")
tmp<-RL1; tmp$beta[,c(3,4,7,8)]<-0
survFRSDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRSDelta)[1:2]<-c("pred","survived")
tmp<-RL1; tmp$beta[,c(1,3,4,5,7,8)]<-0
survDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDelta)[1:2]<-c("pred","survived")

survFrame%<>%rbind(cbind(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL1",300),
                                    Plot=rep("All $\\beta$ Linear Predictor Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDemog,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL1",300),
                                    Plot=rep("Gompertz Demographic Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRS,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL1",300),
                                    Plot=rep("Systolic Mean Term Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRSDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL1",300),
                                    Plot=rep("Mean and $\\Delta$ Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL1",300),
                                    Plot=rep("$\\Delta$ Terms",300))))

survTot<-GetSurvival(RL=RL2,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")
tmp<-RL2; tmp$beta[]<-0
survDemog<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDemog)[1:2]<-c("pred","survived")
tmp<-RL2; tmp$beta[,2:8]<-0
survFRS<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRS)[1:2]<-c("pred","survived")
tmp<-RL2; tmp$beta[,c(3,4,7,8)]<-0
survFRSDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRSDelta)[1:2]<-c("pred","survived")
tmp<-RL2; tmp$beta[,c(1,3,4,5,7,8)]<-0
survDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDelta)[1:2]<-c("pred","survived")

survFrame%<>%rbind(cbind(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL2",300),
                                    Plot=rep("All $\\beta$ Linear Predictor Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDemog,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL2",300),
                                    Plot=rep("Gompertz Demographic Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRS,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL2",300),
                                    Plot=rep("Systolic Mean Term Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRSDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL2",300),
                                    Plot=rep("Mean and $\\Delta$ Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL2",300),
                                    Plot=rep("$\\Delta$ Terms",300))))

survTot<-GetSurvival(RL=RL3,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")
tmp<-RL3; tmp$beta[]<-0
survDemog<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDemog)[1:2]<-c("pred","survived")
tmp<-RL3; tmp$beta[,2:8]<-0
survFRS<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRS)[1:2]<-c("pred","survived")
tmp<-RL3; tmp$beta[,c(3,4,7,8)]<-0
survFRSDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRSDelta)[1:2]<-c("pred","survived")
tmp<-RL3; tmp$beta[,c(1,3,4,5,7,8)]<-0
survDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDelta)[1:2]<-c("pred","survived")

survFrame%<>%rbind(cbind(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL3",300),
                                    Plot=rep("All $\\beta$ Linear Predictor Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDemog,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL3",300),
                                    Plot=rep("Gompertz Demographic Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRS,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL3",300),
                                    Plot=rep("Systolic Mean Term Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRSDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL3",300),
                                    Plot=rep("Mean and $\\Delta$ Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL3",300),
                                    Plot=rep("$\\Delta$ Terms",300))))

survTot<-GetSurvival(RL=RL4,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")
tmp<-RL4; tmp$beta[]<-0
survDemog<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDemog)[1:2]<-c("pred","survived")
tmp<-RL4; tmp$beta[,2:8]<-0
survFRS<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRS)[1:2]<-c("pred","survived")
tmp<-RL4; tmp$beta[,c(3,4,7,8)]<-0
survFRSDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRSDelta)[1:2]<-c("pred","survived")
tmp<-RL4; tmp$beta[,c(1,3,4,5,7,8)]<-0
survDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDelta)[1:2]<-c("pred","survived")

survFrame%<>%rbind(cbind(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL4",300),
                                    Plot=rep("All $\\beta$ Linear Predictor Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDemog,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL4",300),
                                    Plot=rep("Gompertz Demographic Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRS,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL4",300),
                                    Plot=rep("Systolic Mean Term Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRSDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL4",300),
                                    Plot=rep("Mean and $\\Delta$ Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL4",300),
                                    Plot=rep("$\\Delta$ Terms",300))))

survTot<-GetSurvival(RL=RL5,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")
tmp<-RL5; tmp$beta[]<-0
survDemog<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDemog)[1:2]<-c("pred","survived")
tmp<-RL5; tmp$beta[,2:7]<-0
survFRS<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRS)[1:2]<-c("pred","survived")
tmp<-RL5; tmp$beta[,c(3,4,7)]<-0
survFRSDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRSDelta)[1:2]<-c("pred","survived")
tmp<-RL5; tmp$beta[,c(1,3,4,6,7)]<-0
survDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDelta)[1:2]<-c("pred","survived")

survFrame%<>%rbind(cbind(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL5",300),
                                    Plot=rep("All $\\beta$ Linear Predictor Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDemog,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL5",300),
                                    Plot=rep("Gompertz Demographic Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRS,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL5",300),
                                    Plot=rep("FRS-1998 Term Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRSDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL5",300),
                                    Plot=rep("FRS-1998 and $\\Delta$ Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL5",300),
                                    Plot=rep("$\\Delta$ Terms",300))))

survTot<-GetSurvival(RL=RL6,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")
tmp<-RL6; tmp$beta[]<-0
survDemog<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDemog)[1:2]<-c("pred","survived")
tmp<-RL6; tmp$beta[,2:7]<-0
survFRS<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRS)[1:2]<-c("pred","survived")
tmp<-RL6; tmp$beta[,c(3,4,7)]<-0
survFRSDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRSDelta)[1:2]<-c("pred","survived")
tmp<-RL6; tmp$beta[,c(1,3,4,6,7)]<-0
survDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDelta)[1:2]<-c("pred","survived")

survFrame%<>%rbind(cbind(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL6",300),
                                    Plot=rep("All $\\beta$ Linear Predictor Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDemog,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL6",300),
                                    Plot=rep("Gompertz Demographic Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRS,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL6",300),
                                    Plot=rep("FRS-1998 Term Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRSDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL6",300),
                                    Plot=rep("FRS-1998 and $\\Delta$ Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL6",300),
                                    Plot=rep("$\\Delta$ Terms",300))))

survTot<-GetSurvival(RL=RL7,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")
tmp<-RL7; tmp$beta[]<-0
survDemog<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDemog)[1:2]<-c("pred","survived")
tmp<-RL7; tmp$beta[,2:7]<-0
survFRS<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRS)[1:2]<-c("pred","survived")
tmp<-RL7; tmp$beta[,c(3,4,7)]<-0
survFRSDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRSDelta)[1:2]<-c("pred","survived")
tmp<-RL7; tmp$beta[,c(1,3,4,6,7)]<-0
survDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDelta)[1:2]<-c("pred","survived")

survFrame%<>%rbind(cbind(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL7",300),
                                    Plot=rep("All $\\beta$ Linear Predictor Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDemog,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL7",300),
                                    Plot=rep("Gompertz Demographic Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRS,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL7",300),
                                    Plot=rep("FRS-1998 Term Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRSDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL7",300),
                                    Plot=rep("FRS-1998 and $\\Delta$ Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL7",300),
                                    Plot=rep("$\\Delta$ Terms",300))))

survTot<-GetSurvival(RL=RL8,roc=T,usexhat=F)
names(survTot)[1:2]<-c("pred","survived")
tmp<-RL8; tmp$beta[]<-0
survDemog<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDemog)[1:2]<-c("pred","survived")
tmp<-RL8; tmp$beta[,2:7]<-0
survFRS<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRS)[1:2]<-c("pred","survived")
tmp<-RL8; tmp$beta[,c(3,4,7)]<-0
survFRSDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survFRSDelta)[1:2]<-c("pred","survived")
tmp<-RL8; tmp$beta[,c(1,3,4,6,7)]<-0
survDelta<-GetSurvival(RL=tmp,roc=T,usexhat=F)
names(survDelta)[1:2]<-c("pred","survived")

survFrame%<>%rbind(cbind(calculate_roc(df=survTot,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL8",300),
                                    Plot=rep("All $\\beta$ Linear Predictor Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDemog,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL8",300),
                                    Plot=rep("Gompertz Demographic Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRS,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL8",300),
                                    Plot=rep("FRS-1998 Term Only",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survFRSDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL8",300),
                                    Plot=rep("FRS-1998 and $\\Delta$ Terms",300))))
survFrame%<>%rbind(cbind(calculate_roc(df=survDelta,cost_of_fp =1,cost_of_fn=1,n = 300),
                         data.frame(RL=rep("RL8",300),
                                    Plot=rep("$\\Delta$ Terms",300))))

saveRDS(survFrame,"./Plots/Survival/ROC_Data.Rdata")

survFrame%>%ggplot()+geom_line(aes(fpr,tpr,colour=RL),size=1)+facet_wrap(~Plot,nrow = 2)

survFrame%>%ggplot()+geom_line(aes(fpr,tpr,colour=Plot),size=1)+facet_wrap(~RL,nrow = 2)

survFrame$Plot<-factor(survFrame$Plot)
levels(survFrame$Plot)<-c("Delta Terms","All","FRS and Delta","FRS Only","Gompertz Only","Mean and Delta","Systolic Mean Only")
survFrame$Plot<-as.character(survFrame$Plot)

p<-survFrame%>%ggplot()+geom_point(aes(fpr,tpr,colour=Plot,shape=Event),size=1)+
  scale_color_discrete(labels=TeX(unique(survFrame$Plot))) + geom_abline(slope = 1,intercept = 0) +
  xlab("False Positive Ratio") + ylab("True Positive Ratio") +
  facet_wrap(~RL,labeller = as_labeller(TeX,default = label_parsed),nrow = 2);p
ggsave("ROCSurvival_RL.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 12,height = 6) 

p<-survFrame%>%ggplot()+geom_point(aes(fpr,tpr,colour=RL,shape=Event),size=1)+
  xlab("False Positive Ratio") + ylab("True Positive Ratio") + geom_abline(slope = 1,intercept = 0) +
  # ggtitle(label = c("Delta Terms","All","FRS and Delta","FRS Only","Gompertz Only","Mean and Delta","Systolic Mean Only")) +
  facet_wrap(~Plot,nrow = 2);p
ggsave("ROCSurvival_Plot.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 12,height = 6) 

survFrame$FRS<-F
survFrame$FRS[survFrame$RL%in%c("RL5","RL6","RL7","RL8")]<-T
survFrame$Model<-survFrame$Plot
survFrame$Model[survFrame$Plot%in%c("FRS and Delta","Mean and Delta")]<-"FRS/Mean + Delta"
survFrame$Model[survFrame$Plot%in%c("FRS Only","Systolic Mean Only")]<-"FRS/Sys-Mean Only"

survFrame$Event[survFrame$Event=="All"]<-"All Deaths"
survFrame$Event[survFrame$Event=="CVDHrt"]<-"CVD & Hrt Only"

AUROC<-survFrame%>%filter(RL%in%c("RL1","RL2","RL5","RL6"))%>%group_by(RL,Model,Event,FRS)%>%summarise(AUROC=max(auroc),.groups = "keep")
AUROC$label<-paste0("AUC=",round(AUROC$AUROC,2))
AUROC$x<-0.65
AUROC$y<-0.2
AUROC$y[AUROC$FRS]<-0.3

p<-survFrame%>%filter(RL%in%c("RL1","RL2","RL5","RL6"))%>%ggplot()+geom_line(aes(fpr,tpr,colour=FRS,linetype=FRS),size=1)+
  geom_abline(slope = 1,intercept = 0) + geom_text(data=AUROC,aes(x,y,label=label,colour=FRS))+
  xlab("False Positive Ratio") + ylab("True Positive Ratio") +
  facet_wrap(Event~Model,nrow = 2);p
ggsave("ROCSurvival_Model.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 12,height = 6) 

# He wanted plots with FRS-population for FRS against Mean, then including other terms
# Should have one plot that compares each of them together: FRS vs Mean, all(FRS) vs all(Mean), Gompertz(FRS) vs Gompertz(FRS), Delta(FRS) vs Delta(Mean)
# Shape should be whether the model uses FRS or not
# Could do 2 plots, one for CVDHrt and one for all


output<-data.frame()
for(year in ceiling(min(list_nhanesFRS$T)):ceiling(max(list_nhanesFRS$T)+2)){
  # for(age in ceiling(min(list_nhanesFRS$T+list_nhanesFRS$age)):ceiling(max(list_nhanesFRS$T+list_nhanesFRS$age)+2)){
  for(age in 50:99){  
     
    indies<-list_nhanesFRS$T<=(year+1) & list_nhanesFRS$T>year &
            (list_nhanesFRS$T+list_nhanesFRS$age)<=(age+1) & (list_nhanesFRS$T+list_nhanesFRS$age)>age
    
    output%<>%rbind(data.frame(age=(age+1),
                                sheep.yr=(year+1),
                                total=max((calculate_roc(df=survTot[indies,],cost_of_fp =1,cost_of_fn=1,n = 300))$auroc),
                                Demog=max((calculate_roc(df=survDemog[indies,],cost_of_fp =1,cost_of_fn=1,n = 300))$auroc),
                                FRSonly=max((calculate_roc(df=survFRS[indies,],cost_of_fp =1,cost_of_fn=1,n = 300))$auroc),
                                FRSDelta=max((calculate_roc(df=survFRSDelta[indies,],cost_of_fp =1,cost_of_fn=1,n = 300))$auroc),
                                Delta=max((calculate_roc(df=survDelta[indies,],cost_of_fp =1,cost_of_fn=1,n = 300))$auroc),
                                count=sum(indies),
                                deaths=sum(list_nhanesFRS$eventCVDHrt[indies])))
  }
}

total<-output[!is.na(output$total),]%>%PrepLexisData
tmp<-output[,c(1,2,8)]; names(tmp)[3]<-"total"
countys<-tmp[!is.na(tmp$total),]%>%PrepLexisData
tmp<-output; names(tmp)[c(3,9)]<-c("tmp","total")
deaths<-tmp[!is.na(tmp$total),]%>%PrepLexisData


countys<-tmp[!is.na(tmp$total),]%>%PrepLexisData(funcy=function(i) log(i+1))
tmp<-output; names(tmp)[c(3,9)]<-c("tmp","total")
deaths<-tmp[!is.na(tmp$total),]%>%PrepLexisData(funcy=function(i) log(i+1))


tmp<-output; names(tmp)[c(3,4)]<-c("tmp","total")
demog<-tmp[!is.na(tmp$total),]%>%PrepLexisData
tmp<-output; names(tmp)[c(3,5)]<-c("tmp","total")
FRSonly<-tmp[!is.na(tmp$total),]%>%PrepLexisData
tmp<-output; names(tmp)[c(3,6)]<-c("tmp","total")
FRSDelta<-tmp[!is.na(tmp$total),]%>%PrepLexisData
tmp<-output; names(tmp)[c(3,7)]<-c("tmp","total")
Delta<-tmp[!is.na(tmp$total),]%>%PrepLexisData

p<-lexis_polygon(lg = mylexis, x = total$date, y = total$age, group = total$id,fill = total$cols) +
  labs(x="T - Time Since Census Start",y="Age at Time of Outcome Measured",title = "Area Under ROC") +theme(plot.title = element_text(hjust = 0.5));p
ggsave("LexisAll.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 
cols<-RColorBrewer::brewer.pal(n = 11, name = "Spectral") ; cols<-colorRampPalette(cols)(19)
q<-plot_discrete_cbar(round(unique((histss(na.omit(output$total),n = 20,plotting = F))$breaks),2), 
                      spacing_scaling =3, spacing = "constant", colors = cols,legend_direction = "vertical",direction = 1)
ggsave("LexisAll_col.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 

p<-lexis_polygon(lg = mylexis, x = countys$date, y = countys$age, group = countys$id,fill = countys$cols) +
  labs(x="T - Time Since Census Start",y="Age at Time of Outcome Measured",title = "Number of Individuals") +theme(plot.title = element_text(hjust = 0.5));p
ggsave("LexisCounts.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 
# p<-plot_discrete_cbar(round(unique((histss(output$count,10,plot = F))$breaks),1), 
#                       spacing_scaling = 3,spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = 1)
cols<-RColorBrewer::brewer.pal(n = 11, name = "Spectral") ; cols<-colorRampPalette(cols)(19)
p<-plot_discrete_cbar(round(unique((histss(na.omit(output$count),n = 20,plotting = F))$breaks),2), 
                      spacing_scaling =3, spacing = "constant", colors = cols,legend_direction = "vertical",direction = 1)
ggsave("LexisCounts_col.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 

p<-lexis_polygon(lg = mylexis, x = deaths$date, y = deaths$age, group = deaths$id,fill = deaths$cols) +
  labs(x="T - Time Since Census Start",y="Age at Time of Outcome Measured",title = "Number of Deaths") +theme(plot.title = element_text(hjust = 0.5));p
ggsave("LexisDeaths.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 
# p<-plot_discrete_cbar(round(unique((histss(output$deaths,10,plot = F))$breaks),1), 
#                       spacing_scaling = 3,spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = 1)
cols<-RColorBrewer::brewer.pal(n = 11, name = "Spectral") ; cols<-colorRampPalette(cols)(19)
p<-plot_discrete_cbar(round(unique((histss(na.omit(output$deaths),n = 20,plotting = F))$breaks),2), 
                      spacing_scaling =3, spacing = "constant", colors = cols,legend_direction = "vertical",direction = 1)
ggsave("LexisDeaths_col.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 




p<-lexis_polygon(lg = mylexis, x = demog$date, y = demog$age, group = demog$id,fill = demog$cols) +
  labs(x="T - Time Since Census Start",y="Age at Time of Outcome Measured",title = "") +theme(plot.title = element_text(hjust = 0.5));p
ggsave("LexisDemog.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 
p<-plot_discrete_cbar(unique((hist(output$demog,breaks = 10,plot = F))$breaks), 
                      spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = -1)
ggsave("LexisDemog_col.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 

p<-lexis_polygon(lg = mylexis, x = FRSonly$date, y = FRSonly$age, group = FRSonly$id,fill = FRSonly$cols) +
  labs(x="T - Time Since Census Start",y="Age at Time of Outcome Measured",title = "") +theme(plot.title = element_text(hjust = 0.5));p
ggsave("LexisFRS.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 
p<-plot_discrete_cbar(unique((hist(output$FRSonly,breaks = 10,plot = F))$breaks), 
                      spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = -1)
ggsave("LexisFRS_col.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 

p<-lexis_polygon(lg = mylexis, x = FRSDelta$date, y = FRSDelta$age, group = FRSDelta$id,fill = FRSDelta$cols) +
  labs(x="T - Time Since Census Start",y="Age at Time of Outcome Measured",title = "") +theme(plot.title = element_text(hjust = 0.5));p
ggsave("LexisFRSDelta.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 
p<-plot_discrete_cbar(unique((hist(output$FRSDelta,breaks = 10,plot = F))$breaks), 
                      spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = -1)
ggsave("LexisFRSDelta_col.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 

p<-lexis_polygon(lg = mylexis, x = Delta$date, y = Delta$age, group = Delta$id,fill = Delta$cols) +
  labs(x="T - Time Since Census Start",y="Age at Time of Outcome Measured",title = "") +theme(plot.title = element_text(hjust = 0.5));p
ggsave("LexisDelta.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 
p<-plot_discrete_cbar(unique((hist(output$Delta,breaks = 10,plot = F))$breaks), 
                      spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = -1)
ggsave("LexisDelta_col.png", plot=p,path = paste0(directory,'Plots/Survival'),width = 6,height = 10) 

# BEST SURVIVAL PREDICTION
# survy<-GetSurvival(RL5)
# j<-which.min(rowMeans(abs(survy$H_j-survy$d_j)))
# RL5$beta[j,]
# colMeans(RL5$gompertz)
# ggplot(data.frame(H_j=survy$H_j[j,],d_j=survy$d_j))+geom_point(aes(H_j,d_j))+geom_abline(slope = 1)

# Fxhat_A<-function(xhat,beta,M_i_S,M_i_D,D_i_S,D_i_D,tau_C_S,tau_C_D,tau_H_S,tau_H_D,nod=F){
#   
#   dimmy<-dim(D_i_S)
#   
#   # beta[ 1:1000 ]
#   # x[ 1:1000  , 1:18018 ]
#   A1<-beta[,1]*(M_i_S - xhat[1,1])
#   A2<-beta[,2]*(abs(D_i_S) - xhat[2,1])
#   A3<-beta[,3]*(tau_C_S - xhat[3,1])
#   A4<-beta[,4]*(tau_H_S - xhat[4,1])
#   A5<-beta[,5]*(M_i_D - xhat[1,2])
#   A6<-beta[,6]*(abs(D_i_D) - xhat[2,2])
#   A7<-beta[,7]*(tau_C_D - xhat[3,2])
#   A8<-beta[,8]*(tau_H_D - xhat[4,2])  
#   
#   B<-exp( A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 )
#   
#   # Save RAM by replacing A values with exp(A)
#   A1<-exp(A1)
#   A2<-exp(A2)
#   A3<-exp(A3)
#   A4<-exp(A4)
#   A5<-exp(A5)
#   A6<-exp(A6)
#   A7<-exp(A7)
#   A8<-exp(A8)
#   
#   C<-rowSums(B)-dimmy[2]
#   
#   C1<-rowSums(A1)-dimmy[2]
#   C2<-rowSums(A2)-dimmy[2]
#   C3<-rowSums(A3)-dimmy[2]
#   C4<-rowSums(A4)-dimmy[2]
#   C5<-rowSums(A5)-dimmy[2]
#   C6<-rowSums(A6)-dimmy[2]
#   C7<-rowSums(A7)-dimmy[2]
#   C8<-rowSums(A8)-dimmy[2]
#   
#   D<-mean(abs(C))
#   
#   D1<-mean(abs(C1))
#   D2<-mean(abs(C2))
#   D3<-mean(abs(C3))
#   D4<-mean(abs(C4))
#   D5<-mean(abs(C5))
#   D6<-mean(abs(C6))
#   D7<-mean(abs(C7))
#   D8<-mean(abs(C8))
#   
#   funcy<-8*D+D1+D2+D3+D4+D5+D6+D7+D8
#     
#   if(nod) return(funcy)
#   
#   derivy<-rep(0,8)
#   
#   for (k in 1:8){
#     
#     derivy[k]<-8*mean(-beta[,k]*(C+dimmy[2])*sign(C)) + mean(-beta[,k]*(get(paste0("C",k))+dimmy[2])*sign(get(paste0("C",k))))
#     # derivy[k]<-( C/D*8*mean( -beta[,k]*B ) + get(paste0("C",k))/get(paste0("D",k))*mean(-beta[,k]*get(paste0("A",k))))
#     
#   }
#   return(array(funcy/derivy,dim = c(4,2)))
# 
# }
# 
# Fxhat_FRS<-function(xhat,beta,FRS,D_i_S,D_i_D,tau_C_S,tau_C_D,tau_H_S,tau_H_D,nod=F){
#   
#   dimmy<-dim(D_i_S)
#   
#   # beta[ 1:1000 ]
#   # x[ 1:1000  , 1:18018 ]
#   A1<-beta[,1]%o%FRS
#   A2<-beta[,2]*(abs(D_i_S) - xhat[2,1])
#   A3<-beta[,3]*(tau_C_S - xhat[3,1])
#   A4<-beta[,4]*(tau_H_S - xhat[4,1])
#   A5<-beta[,5]*(abs(D_i_D) - xhat[2,2])
#   A6<-beta[,6]*(tau_C_D - xhat[3,2])
#   A7<-beta[,7]*(tau_H_D - xhat[4,2])
#   
#   B<-exp( A1 + A2 + A3 + A4 + A5 + A6 + A7)
#   
#   # Save RAM by replacing A values with exp(A)
#   rm(A1)
#   A2<-exp(A2)
#   A3<-exp(A3)
#   A4<-exp(A4)
#   A5<-exp(A5)
#   A6<-exp(A6)
#   A7<-exp(A7)
#   
#   C<-rowSums(B)-dimmy[2]
#   
#   C2<-rowSums(A2)-dimmy[2]
#   C3<-rowSums(A3)-dimmy[2]
#   C4<-rowSums(A4)-dimmy[2]
#   C5<-rowSums(A5)-dimmy[2]
#   C6<-rowSums(A6)-dimmy[2]
#   C7<-rowSums(A7)-dimmy[2]
#   
#   D<-mean(abs(C))
#   
#   D2<-mean(abs(C2))
#   D3<-mean(abs(C3))
#   D4<-mean(abs(C4))
#   D5<-mean(abs(C5))
#   D6<-mean(abs(C6))
#   D7<-mean(abs(C7))
#   
#   funcy<-6*D+D2+D3+D4+D5+D6+D7
#   
#   if(nod) return(funcy)
#   
#   derivy<-rep(0,6)
#   
#   for (k in 2:7){
#     
#     derivy[k-1]<-6*mean(-beta[,k]*(C+dimmy[2])*sign(C)) + mean(-beta[,k]*(get(paste0("C",k))+dimmy[2])*sign(get(paste0("C",k))))
#     # derivy[k-1]<-( C/D*6*mean( -beta[,k]*B ) + get(paste0("C",k))/get(paste0("D",k))*mean(-beta[,k]*get(paste0("A",k))))
#     
#   }
#   
#   tmp<-array(rep(0,8),dim = c(4,2))
#   tmp[2:4,1]<-funcy/derivy[1:3]
#   tmp[2:4,2]<-funcy/derivy[4:6]
#   
#   return(tmp)
#   
# }
