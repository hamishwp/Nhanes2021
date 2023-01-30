library(magrittr)
# Define the model
FRSt<-F
FRS1998<-T
CVDHrt<-T
if(FRSt){
  # Extract the appropriate data
  lister<-list_nhanesFRS
  nhanes<-DF_Nhanes(list_nhanesFRS)
  # Setup the initial values
  if(CVDHrt) {
    initials<-c(log(c(apply(RL7$gompertz,2:4,mean))),colMeans(RL7$beta))
  } else {
    initials<-c(log(c(apply(RL8$gompertz,2:4,mean))),colMeans(RL8$beta))
  }
} else {
  # Extract the appropriate data
  lister<-list_nhanesA
  nhanes<-DF_Nhanes(list_nhanesA)
  # Setup the initial values
  if(CVDHrt) {
    initials<-c(log(c(apply(RL1$gompertz,2:4,mean))),colMeans(RL1$beta))
  } else {
    initials<-c(log(c(apply(RL2$gompertz,2:4,mean))),colMeans(RL2$beta))
  }
}
# Extract the individual-level covariates:
M_i_S<-rowMeans(cbind(as.matrix(lister$sys),as.matrix(lister$sys_home)))
M_i_D<-rowMeans(cbind(as.matrix(lister$dias),as.matrix(lister$dias_home)))
D_i_S<-abs(rowMeans(lister$sys)-rowMeans(lister$sys_home))
D_i_D<-abs(rowMeans(lister$dias)-rowMeans(lister$dias_home))
tau_C_S<-apply(lister$sys,1,sd)
tau_C_D<-apply(lister$dias,1,sd)
tau_H_S<-apply(lister$sys_home,1,sd)
tau_H_D<-apply(lister$dias_home,1,sd)
age<-nhanes$age
Time<-nhanes$T  
if(CVDHrt) delta<-nhanes$eventCVDHrt else delta<-nhanes$eventall
gender<-nhanes$male+1
ethn<-nhanes$black+2*nhanes$white+3*nhanes$other
LLL<-length(D_i_S)
# Auto-centering function
fxhat<-function(x,beta) (log(sum(exp(x*beta))) - log(length(x)) ) / beta
# Linear predictor terms depend on the model chosen 
if(FRSt){
  if(FRS1998) FRS<-lister$FRS.1998 else FRS<-lister$FRS.ATP
  # If using the Framingham Risk Score
  linpred<-function(beta){
    A1<-beta[1]*(FRS-fxhat(FRS,beta[1]))
    A2<-beta[2]*(abs(D_i_S)  - fxhat(abs(D_i_S),beta[2]))
    A3<-beta[3]*(tau_C_S     - fxhat(tau_C_S,beta[3]))
    A4<-beta[4]*(tau_H_S     - fxhat(tau_H_S,beta[4]))
    A5<-beta[5]*(abs(D_i_D)  - fxhat(abs(D_i_D),beta[5]))
    A6<-beta[6]*(tau_C_D     - fxhat(tau_C_D,beta[6]))
    A7<-beta[7]*(tau_H_D     - fxhat(tau_H_D,beta[7])) 
    # total linear predictor
    A1 + A2 + A3 + A4 + A5 + A6 + A7
  }
}   else {
  # Otherwise
  linpred<-function(beta){
    A1<-beta[1]*(M_i_S       - fxhat(M_i_S,beta[1]))
    A2<-beta[2]*(abs(D_i_S)  - fxhat(abs(D_i_S),beta[2]))
    A3<-beta[3]*(tau_C_S     - fxhat(tau_C_S,beta[3]))
    A4<-beta[4]*(tau_H_S     - fxhat(tau_H_S,beta[4]))
    A5<-beta[5]*(M_i_D       - fxhat(M_i_D,beta[5]))
    A6<-beta[6]*(abs(D_i_D)  - fxhat(abs(D_i_D),beta[6]))
    A7<-beta[7]*(tau_C_D     - fxhat(tau_C_D,beta[7]))
    A8<-beta[8]*(tau_H_D     - fxhat(tau_H_D,beta[8]))  
    # total linear predictor
    A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8
  }
}  
# Likelihood function
Likely<-function(input){
  # Convert from input parameter space value to covariates
  gompertz<-array(exp(input[1:12]), c(2,3,2))
  beta<-input[13:length(input)]
  # Sort gompertz into B and theta
  Btheta<-t(vapply(1:LLL,FUN = function(j) gompertz[gender[j],ethn[j],],
                 FUN.VALUE = numeric(2)))
  # Calculate the hazard risk
  H_j<-Btheta[,1]*exp(Btheta[,2]*(age + Time) + linpred(beta))
  # H_j<-(Btheta[,1]/Btheta[,2])*exp(age*Btheta[,2] + linpred(beta))*(exp(Btheta[,2]*Time) - 1)
  # Then the global concordance index
  concordance.index(x = H_j,surv.time = age+Time,
                    surv.event = delta)$c.index
}

# Check the results from the likelihood function just in case
Likely(initials)

# Reduce the RAM required during simulations
includers<-c("fxhat","linpred","Likely","initials","LLL",
             "D_i_S","D_i_D","tau_C_D","tau_C_S","tau_H_D","tau_H_S",
             "gender","ethn","age","Time","delta")
if(FRSt) includers%<>%c("FRS") else includers%<>%c("M_i_S","M_i_D")
rm(list=setdiff(ls(), includers))

# Let's go!
output<-optim(fn = Likely,par=initials,
              control = list(fnscale=-1,maxit=5000),hessian = T)      
   
saveRDS(output,"./RL1.RData")

values<-output$par; values[1:12]<-exp(values[1:12])
copinit<-initials; copinit[1:12]<-exp(copinit[1:12])
print("Gompertz:")
values[1:12]
copinit[1:12]
print("Beta:")
values[13:length(values)]
copinit[13:length(values)]
print("C-index:")
output$value

DF<-data.frame(deaths=delta,Time=age+Time,Year=Time)
DF$deaths[DF$deaths==0]<-"Censored"
DF$deaths[DF$deaths==1]<-"Death"

p<-DF%>%ggplot(aes(Time, group=deaths))+
  geom_density(aes(y = ..count..,fill=deaths),position = "stack") +
  xlab("Time (Age + Census Time)") + ylab("Number of Deaths")
ggsave("1D_Density_Events.png",p)

p<-DF%>%ggplot(aes(Time, Year))+stat_density_2d_filled(bins=50,show.legend = FALSE)+
  facet_wrap(~deaths) + xlab("Time (Age + Census Time)") + ylab("Years Since Census")
ggsave("2D_Density_Events.png",p)
