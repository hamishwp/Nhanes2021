library(tidyverse)
library(latex2exp)
library(ggplot2)
library(dplyr)
library(pracma)
library(magrittr)
library(ggpubr)

# directory<-'~/Documents/Coding/Oxford/BP_NhanesA/'
directory<-"/home/patten/Documents/Coding/Oxford/Nhanes2021/"

load(paste0(directory,"Data_cleaned/nhanes_cleaned_lists.RData"))

FRS<-F

DF_Nhanes<-function(list_nhanes){
  
  DF<-as.data.frame.list(list_nhanes[c("eventall","eventCVDHrt","black","white","other",
                                       "female","male","T","age")])
  if(!is.null(list_nhanes[["FRS.ATP"]]) & !is.null(list_nhanes[["FRS.1998"]]))
    DF%<>%cbind(as.data.frame.list(list_nhanes[c("FRS.ATP","FRS.1998")]))
  
  return(DF)
}
if(FRS) nhanes<-DF_Nhanes(list_nhanesFRS) else nhanes<-DF_Nhanes(list_nhanesA)

costfunction<-function(theta,nhanes,inds,delta){
  
  age<-nhanes$age[inds]
  Time<-nhanes$T[inds]
  
  n<-30
  lenny<-n-2
  breaks<-seq(from=min(age,na.rm = T)-1, to=max(age,na.rm = T)+1,length.out=n)
  
  a_j<-d_j<-H_j<-vector(length = lenny)
  
  totH<-sum(exp(age*theta)*(exp(theta*Time) - 1))
  totd<-sum(delta)
  
  for (i in 2:(n-1)){
    
    # Age histogram a_j
    a_j[i-1]<-breaks[i]
    # indices of people that lived through age a_j
    iii<-age<=breaks[i]
    
    # cumulative hazard of all people that lived through age a_j
    H_j[i-1]<-sum(exp(age[iii]*theta)*(exp(theta*Time[iii]) - 1))/totH
    # cumulative deaths of all people that lived up to age a_j
    d_j[i-1]<-sum(delta[iii])/totd
    
  }
  
  return(abs(H_j[lenny]-d_j[lenny])*mean(abs(H_j-d_j)))
  
}

CalcB<-function(theta,nhanes,inds,delta){
  age<-nhanes$age[inds]
  Time<-nhanes$T[inds]
  
  n<-30
  lenny<-n-2
  breaks<-seq(from=min(age,na.rm = T)-1, to=max(age,na.rm = T)+1,length.out=n)
  
  a_j<-d_j<-H_j<-vector(length = lenny)
  
  totH<-sum(exp(age*theta)*(exp(theta*Time) - 1))
  totd<-sum(delta)
  
  for (i in 2:(n-1)){
    
    # Age histogram a_j
    a_j[i-1]<-breaks[i]
    # indices of people that lived through age a_j
    iii<-age<=breaks[i]
    
    # cumulative hazard of all people that lived through age a_j
    H_j[i-1]<-sum(exp(age[iii]*theta)*(exp(theta*Time[iii]) - 1))
    # cumulative deaths of all people that lived up to age a_j
    d_j[i-1]<-sum(delta[iii])
    
  }
  
  return(mean(d_j*theta-H_j))
  
}

# Black Females
bftheta<-optimize(f=function(theta) costfunction(theta = theta,
                                            nhanes = nhanes,
                                            inds = nhanes$female==1 & nhanes$black==1,
                                            delta=nhanes$eventCVDHrt),
             lower=0, upper=0.25); bftheta$minimum

# Black Males
bmtheta<-optimize(f=function(theta) costfunction(theta = theta,
                                            nhanes = nhanes,
                                            inds = nhanes$male==1 & nhanes$black==1,
                                            delta=nhanes$eventCVDHrt),
             lower=0, upper=0.25); bmtheta$minimum

# Other Females
oftheta<-optimize(f=function(theta) costfunction(theta = theta,
                                            nhanes = nhanes,
                                            inds = nhanes$female==1 & nhanes$other==1,
                                            delta=nhanes$eventCVDHrt),
             lower=0, upper=0.25); oftheta$minimum

# Other Males
omtheta<-optimize(f=function(theta) costfunction(theta = theta,
                                            nhanes = nhanes,
                                            inds = nhanes$male==1 & nhanes$other==1,
                                            delta=nhanes$eventCVDHrt),
             lower=0, upper=0.25); omtheta$minimum

# White Females
wftheta<-optimize(f=function(theta) costfunction(theta = theta,
                                            nhanes = nhanes,
                                            inds = nhanes$female==1 & nhanes$white==1,
                                            delta=nhanes$eventCVDHrt),
            lower=0, upper=0.25); wftheta$minimum

# White Males
wmtheta<-optimize(f=function(theta) costfunction(theta = theta,
                                            nhanes = nhanes,
                                            inds = nhanes$male==1 & nhanes$white==1,
                                            delta=nhanes$eventCVDHrt),
             lower=0, upper=0.25); wmtheta$minimum


# Black Females
bfB<-costfunction(theta = bftheta$minimum,
                 nhanes = nhanes,
                 inds = nhanes$female==1 & nhanes$black==1,
                 delta=nhanes$eventCVDHrt);bfB

# Black Males
bmB<-costfunction(theta = bmtheta$minimum,
                  nhanes = nhanes,
                  inds = nhanes$male==1 & nhanes$black==1,
                  delta=nhanes$eventCVDHrt);bmB

# Other Females
ofB<-costfunction(theta = oftheta$minimum,
                  nhanes = nhanes,
                  inds = nhanes$female==1 & nhanes$other==1,
                  delta=nhanes$eventCVDHrt);ofB

# Other Males
omB<-costfunction(theta = omtheta$minimum,
                  nhanes = nhanes,
                  inds = nhanes$male==1 & nhanes$other==1,
                  delta=nhanes$eventCVDHrt);omB

# White Females
wfB<-costfunction(theta = wftheta$minimum,
                  nhanes = nhanes,
                  inds = nhanes$female==1 & nhanes$white==1,
                  delta=nhanes$eventCVDHrt);wfB

# White Males
wmB<-costfunction(theta = wmtheta$minimum,
                  nhanes = nhanes,
                  inds = nhanes$male==1 & nhanes$white==1,
                  delta=nhanes$eventCVDHrt);wmB

print(paste0("Black Female: Deaths = ",
             signif(sum(nhanes$eventCVDHrt[nhanes$female==1 & nhanes$black==1])
             /sum(nhanes$female==1 & nhanes$black==1),3),
             ", Gompertz B = ",signif(bfB,3),
             ", Gompertz theta = ",signif(bftheta$minimum,3)))

print(paste0("Black Male: Deaths = ",
             signif(sum(nhanes$eventCVDHrt[nhanes$male==1 & nhanes$black==1])
                    /sum(nhanes$male==1 & nhanes$black==1),3),
             ", Gompertz B = ",signif(bmB,3),
             ", Gompertz theta = ",signif(bmtheta$minimum,3)))

print(paste0("Other Female: Deaths = ",
             signif(sum(nhanes$eventCVDHrt[nhanes$female==1 & nhanes$other==1])
                    /sum(nhanes$female==1 & nhanes$other==1),3),
             ", Gompertz B = ",signif(ofB,3),
             ", Gompertz theta = ",signif(oftheta$minimum,3)))

print(paste0("Other Male: Deaths = ",
             signif(sum(nhanes$eventCVDHrt[nhanes$male==1 & nhanes$other==1])
                    /sum(nhanes$male==1 & nhanes$other==1),3),
             ", Gompertz B = ",signif(omB,3),
             ", Gompertz theta = ",signif(omtheta$minimum,3)))

print(paste0("White Female: Deaths = ",
             signif(sum(nhanes$eventCVDHrt[nhanes$female==1 & nhanes$white==1])
                    /sum(nhanes$female==1 & nhanes$white==1),3),
             ", Gompertz B = ",signif(wfB,3),
             ", Gompertz theta = ",signif(wftheta$minimum,3)))

print(paste0("White Male: Deaths = ",
             signif(sum(nhanes$eventCVDHrt[nhanes$male==1 & nhanes$white==1])
                    /sum(nhanes$male==1 & nhanes$white==1),3),
             ", Gompertz B = ",signif(wmB,3),
             ", Gompertz theta = ",signif(wmtheta$minimum,3)))


# # Old NHANES
# # bftheta$minimum
# 0.05411574
# # bmtheta$minimum
# 0.008821352
# # oftheta$minimum
# 0.01645104
# # omtheta$minimum
# 0.01521158
# # wftheta$minimum
# 0.03604716
# # wmtheta$minimum
# 0.04981642
# 
# 
# # bfB
# 1.930754e-07
# # bmB
# 1.927929e-08
# # ofB
# 1.026715e-07
# # omB
# 5.224348e-09
# # wfB
# 2.783449e-07
# # wmB
# 4.785367e-07
# 
# 
# 
# # New NHANES 2021:
# # bftheta$minimum
# 0.0309959
# # bmtheta$minimum
# 0.01799212
# # oftheta$minimum
# 0.01177535
# # omtheta$minimum
# 0.008114348
# # wftheta$minimum
# 0.03079643
# # wmtheta$minimum
# 0.01679283
# 
# 
# # bfB
# 2.24593e-07
# # bmB
# 1.572613e-08
# # ofB
# 1.227338e-08
# # omB
# 5.071319e-06
# # wfB
# 3.115477e-07
# # wmB
# 9.763834e-08
