functions {
  
  /*
  * Declarations
  */
  real log_h(real linpred, real T, real B, real theta, real age);
  real H_t(real linpred, real T, real B, real theta, real age);
  
  real surv_dens_lpdf(vector T_and_delta, real linpred, real B, real theta, real age){
    
    real T;
    real delta;

    T     = T_and_delta[1];
    delta = T_and_delta[2];
    
    return (delta * log_h(linpred, T, B, theta, age))
    -    H_t(linpred, T, B, theta, age) ;
    
  }
  
  real log_h(real linpred, real T, real B, real theta, real age){
    
    return log(B) + theta*(age + T) + linpred;
    
  }
  real H_t(real linpred, real T, real B, real theta, real age){
    
    return (B/theta)*exp(age*theta + linpred)*(exp(theta*T) - 1);
    
  }

  vector linear_predictor(vector FRS, vector m1, matrix Delta, vector m2,
                          matrix vterm, vector mean_vterm, vector beta){
    
     return((FRS - m1[1]) * beta[1]  +
     (Delta[,1] - m2[1])  * beta[2]  +
     (vterm[,1] - mean_vterm[1]) * beta[3] +
     (vterm[,2] - mean_vterm[2]) * beta[4] +
     (Delta[,2] - m2[2]) * beta[5]  +
     (vterm[,3] - mean_vterm[3]) * beta[6] +
     (vterm[,4] - mean_vterm[4]) * beta[7]);
     
  }
  
}

data {
  
  int<lower=0> N;
  int<lower=0> max_nj;
  
  real yhomeS[N,max_nj];
  real yclinicS[N,max_nj];
  real yhomeD[N,max_nj];
  real yclinicD[N,max_nj];
  
  vector[N] FRS;
  vector[N] age;
  int<lower=1,upper=2>  sex[N];
  int<lower=1,upper=3> race[N];
  
  vector[2] Surv_T_and_delta[N];
  
  // The EmpBayes Parameters
  vector[2] m_M;
  vector[2] tau_Delta;
  vector[2] tau_M;
  vector[2] m_Delta;
  
  vector[2] alpha_HOME;
  vector[2] alpha_CLINIC;
  vector[2] beta_HOME;
  vector[2] beta_CLINIC;

  vector[2] FRSc;

  vector[2] xhat_M;
  vector[2] xhat_D;
  vector[4] xhat_v;

}
parameters {

 /*
 *  See section 22.5 in the Stan manual for an explanation of array and vector referencing.
 *  Basically:
 *  vector[i] thing[j,k,l]
 *  is referenced thing[j,k,l,i]
 *  That is, thing[j,k,l] is a vector[i]
 */

  matrix[N,2] M_i; 
  matrix[N,2] Delta_i;

  matrix<lower=0>[N,4] tauis;

  vector[7] beta;
  real<lower=0> gompertz[2,3,2]; // gompertz[sex, race, ]

}
transformed parameters {
  // These are the parameters Stan needs

  matrix<lower=0>[N,2] abs_Delta_i;
  vector[N] B;
  vector[N] theta;
  vector[N] linnie;
  // vector[2] xhat_M;
  // vector[2] xhat_D;
  // vector[4] xhat_v;

  matrix<lower=0>[N,4] sigmais;

  sigmais = 1. ./ sqrt(tauis);
  
  abs_Delta_i = fabs(Delta_i);
  
  for(i in 1:N){
    B[i]=gompertz[sex[i], race[i],1];
    theta[i]=gompertz[sex[i], race[i],2];
  }
  
  // xhat_D[1] =  (log_sum_exp(Delta_i[,1]*beta[2]) - log(N)) ./beta[2];
  // xhat_D[2] =  (log_sum_exp(Delta_i[,2]*beta[5]) - log(N)) ./beta[5];

  // xhat_v[1] =  (log_sum_exp(sigmais[,1]*beta[3]) - log(N)) ./beta[3];
  // xhat_v[2] =  (log_sum_exp(sigmais[,2]*beta[4]) - log(N)) ./beta[4];
  // xhat_v[3] =  (log_sum_exp(sigmais[,3]*beta[6]) - log(N)) ./beta[6];
  // xhat_v[4] =  (log_sum_exp(sigmais[,4]*beta[7]) - log(N)) ./beta[7];
  
  linnie=linear_predictor(FRS, FRSc, abs_Delta_i, xhat_D, sigmais, xhat_v, beta);

}

model {
  
  for(i_sex in 1:2){
    for(i_race in 1:3){
      gompertz[i_sex,i_race,1] ~ normal(0, 2); // B
      gompertz[i_sex,i_race,2] ~ normal(0, 2); // theta
    }
  }
  
  beta ~ cauchy(0,100);

  for(j in 1:max_nj) {
    yclinicS[,j] ~ normal(M_i[,1] - Delta_i[,1], sigmais[,1]);
    yhomeS[,j]   ~ normal(M_i[,1] + Delta_i[,1], sigmais[,2]);
    yclinicD[,j] ~ normal(M_i[,2] - Delta_i[,2], sigmais[,3]);
    yhomeD[,j]   ~ normal(M_i[,2] + Delta_i[,2], sigmais[,4]);
  }

  M_i[,1] ~ normal( m_M[1], 1/sqrt(tau_M[1]) );
  Delta_i[,1] ~ normal( m_Delta[1], 1/sqrt(tau_Delta[1]) );
  tauis[,1] ~ gamma(alpha_CLINIC[1], beta_CLINIC[1]);
  tauis[,2] ~ gamma(alpha_HOME[1], beta_HOME[1]);
  
  M_i[,2] ~ normal( m_M[2], 1/sqrt(tau_M[2]) );
  Delta_i[,2] ~ normal( m_Delta[2], 1/sqrt(tau_Delta[2]) );
  tauis[,3] ~ gamma(alpha_CLINIC[2], beta_CLINIC[2]);
  tauis[,4] ~ gamma(alpha_HOME[2], beta_HOME[2]);
  
  for(i in 1:N){
    Surv_T_and_delta[i] ~ surv_dens(linnie[i], B[i], theta[i], age[i]);
  }
  
}


