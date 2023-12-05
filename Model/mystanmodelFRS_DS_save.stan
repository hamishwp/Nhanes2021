functions {

/*
* Declarations
*/

// Have to declare the functions before I can use them


real log_h(real T,
         real FRS, vector m1, vector Delta, vector m2, vector tauis, vector mean_tauis,
         vector beta, real[] gompertz, real age);

real H_t(real T,
         real FRS, vector m1, vector Delta, vector m2, vector tauis, vector mean_tauis,
         vector beta, real[] gompertz, real age);

/*
* Definitions
*/
real surv_dens_lpdf(vector T_and_delta,
                    real FRS , vector m1, vector Delta, vector m2,
                    vector tauis, vector mean_tauis,
                    vector beta, real[] gompertz, real age){
    real T;
    real delta;

    T     = T_and_delta[1];
    delta = T_and_delta[2];

    return (delta * log_h(T, FRS, m1, Delta, m2, tauis, mean_tauis, beta,  gompertz, age))
                 -    H_t(T, FRS, m1, Delta, m2, tauis, mean_tauis, beta,  gompertz, age) ;

    }

real log_h(real T,
         real FRS, vector m1, vector Delta, vector m2,
         vector tauis, vector mean_tauis,
         vector beta, real[] gompertz, real age){

         real B;
         real theta;

         B     = gompertz[1];
         theta = gompertz[2];


      return log(B) + theta*(age + T) +
                (FRS - m1[1]) * beta[1] +      
                (Delta[1] - m2[1])  * beta[2]  +
                (tauis[1] - mean_tauis[1]) * beta[3] +
                (tauis[2] - mean_tauis[2]) * beta[4] +                
                (Delta[2] - m2[2])  * beta[5]  +
                (tauis[3] - mean_tauis[3]) * beta[6] +
                (tauis[4] - mean_tauis[4]) * beta[7];

}
real H_t(real T,
         real FRS, vector m1, vector Delta, vector m2,
         vector tauis, vector mean_tauis,
         vector beta, real[] gompertz, real age){


         real B;
         real theta;

         B     = gompertz[1];
         theta = gompertz[2];

      return (B/theta)*
              exp( age*theta +
                (FRS - m1[1]) * beta[1]  +      
                (Delta[1] - m2[1])  * beta[2]  +
                (tauis[1] - mean_tauis[1]) * beta[3] +
                (tauis[2] - mean_tauis[2]) * beta[4] +                
                (Delta[2] - m2[2])  * beta[5]  +
                (tauis[3] - mean_tauis[3]) * beta[6] +
                (tauis[4] - mean_tauis[4]) * beta[7]    
                )*(exp(theta*T) - 1);
}

}

data {

int<lower=0> N;
int<lower=0> max_nj;

real yhomeS[N,max_nj];
real yclinicS[N,max_nj];
real yhomeD[N,max_nj];
real yclinicD[N,max_nj];

real FRS[N];
real age[N];
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

vector[2] xhat_M;
vector[2] xhat_D;
vector[4] xhat_t;

}
parameters {

 /*
 *  See section 22.5 in the Stan manual for an explanation of array and vector referencing.
 *  Basically:
 *  vector[i] thing[j,k,l]
 *  is referenced thing[j,k,l,i]
 *  That is, thing[j,k,l] is a vector[i]
 */

  vector[2] M_i[N]; 
  vector[2] Delta_i[N];

  vector<lower=0>[4] tauis[N];

  vector[7] beta;
  real<lower=0> gompertz[2,3,2]; // gompertz[sex, race, ]

}
transformed parameters {
  // These are the parameters Stan needs
  vector<lower=0>[4] sigmais[N];

  for(i in 1:N){
    sigmais[i,1] = 1/sqrt(tauis[i,1]);
    sigmais[i,2] = 1/sqrt(tauis[i,2]);
    sigmais[i,3] = 1/sqrt(tauis[i,3]);
    sigmais[i,4] = 1/sqrt(tauis[i,4]);    
  }

}

model {

  vector[2] abs_Delta_i[N];

  for(i_sex in 1:2){
    for(i_race in 1:3){
      gompertz[i_sex,i_race,1] ~ normal(0, 2); // B
      gompertz[i_sex,i_race,2] ~ normal(0, 2); // theta
    }
  }

  beta ~ cauchy(0,100);

  for(i in 1:N){

    for(j in 1:max_nj) {
      yclinicS[i,j] ~ normal(M_i[i,1] - Delta_i[i,1], sigmais[i,1] );
      yhomeS[i,j]   ~ normal(M_i[i,1] + Delta_i[i,1], sigmais[i,2] );
      yclinicD[i,j] ~ normal(M_i[i,2] - Delta_i[i,2], sigmais[i,3] );
      yhomeD[i,j]   ~ normal(M_i[i,2] + Delta_i[i,2], sigmais[i,4] );      
    }

    M_i[i,1] ~ normal( m_M[1], 1/sqrt(tau_M[1]) );
    Delta_i[i,1] ~ normal( m_Delta[1], 1/sqrt(tau_Delta[1]) );
    tauis[i,1] ~ gamma(alpha_CLINIC[1], beta_CLINIC[1]);
    tauis[i,2] ~ gamma(alpha_HOME[1], beta_HOME[1]);
    abs_Delta_i[i,1] = fabs(Delta_i[i,1]);
    
    M_i[i,2] ~ normal( m_M[2], 1/sqrt(tau_M[2]) );
    Delta_i[i,2] ~ normal( m_Delta[2], 1/sqrt(tau_Delta[2]) );
    tauis[i,3] ~ gamma(alpha_CLINIC[2], beta_CLINIC[2]);
    tauis[i,4] ~ gamma(alpha_HOME[2], beta_HOME[2]);
    abs_Delta_i[i,2] = fabs(Delta_i[i,2]);    
    
    Surv_T_and_delta[i] ~ surv_dens(FRS[i], xhat_M,
                                    abs_Delta_i[i], xhat_D,
                                    tauis[i], xhat_t,
                                    beta,
                                    gompertz[sex[i], race[i]],
                                    age[i]);
  }
}


