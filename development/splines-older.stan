
data {
    
    int<lower=0> p;                    // N. covariates
    int<lower=0> K;                    // N. basis expansions
    int<lower=0> n;                    // N. observations
    // int<lower=0> ntest;                // N. observations (ppc)
    matrix[n, K] x[p];                 // Predictor matrix (spline expansion)
    // matrix[ntest, K] xtest;           // Predictor matrix (ppc)
    vector[n] y;                       // log size at time t+1 
    
    int<lower=0> nplot;                // number of plots
    int<lower=1> plot[n];              // index for plot
    
    int<lower=0> ntree;                // number of trees
    int<lower=1> tree[n];              // index for trees
    // int<lower=1> treetest[ntest];     //index for trees (ppc)
    int<lower=1> plotfortree[ntree];   // plot index for each tree
}

parameters {
    
    real beta0;                          // intercept means
    vector[K] beta[p];                   // other coeff mean
    
    real beta0_p_tilde[nplot];           // plot-level intercepts
    real<lower=0> s_beta0_p;             // plot-level intercept variance
    real beta0_t_tilde[ntree];           // tree-level intercepts
    real<lower=0> s_beta0_t;             // tree-level intercept variance
    real<lower=0> sigma_beta[K];         // spline coefficient variances
    real<lower=0> sigma_y;               // Residual for growth model
    
}

transformed parameters {
    real beta0_p[nplot];                 // plot-level intercepts
    real beta0_t[ntree];                 // tree-level intercepts
    
    for(i in 1:nplot){
        beta0_p[i] = beta0 + s_beta0_p * beta0_p_tilde[i];
    }
    
    for(t in 1:ntree){
        beta0_t[t] = beta0_p[plotfortree[t]] + s_beta0_t * beta0_t_tilde[t];
    }
    
}

model {
    vector[n] mu;
    
    beta0 ~ normal(0, 100);
    beta0_p_tilde ~ normal(0,1);
    beta0_t_tilde ~ normal(0,1);

    for (i in 1:p) {
        beta[i][1] ~ normal(0, 1);
      for (k in 2:K) {
          beta[i][k] ~ normal(beta[i][k-1], sigma_beta[i]);
      }
    }

    s_beta0_p ~ cauchy(0, 2.5);
    s_beta0_t ~ cauchy(0, 2.5);
    
    sigma_beta ~ cauchy(0, 2.5);
    sigma_y ~ gamma(2, 0.01);
    
    //tried nesting random effects of trees within plots--had issue
    //for(p in 1:nplot){
        //beta0_p[p] ~ normal(u_beta0, s_beta0_p);
        //for(t in 1:ntree){
            //beta0_t[t] ~ normal(beta0_p[plotfortree[t]], s_beta0_t);
            //}
            //}
            
            
            // GROWTH MODEL
            
            for(i in 1:n){
                real tmp = 0;
                for (j in 1:p) {
                    tmp = tmp + x[j][i] * beta[j];
                }
                mu[i] = beta0_t[tree[i]] + tmp;
            }
            
            y ~ normal(mu, sigma_y);
            //yG ~ gamma(mG,sigma_y);
            
}

// generated quantities{
//     vector[nG] yrep;
//     vector[nGtest] ypred;
//     vector[nG] log_lik;
//     
//     for(n in 1:nG){
//         yrep[n] = normal_rng(beta0_t[tree[n]] + xG[n] * u_beta, sigma_y);
//         log_lik[n] = normal_lpdf(yG[n] | beta0_t[tree[n]] + xG[n] * u_beta, sigma_y);
//     }
//     for(n in 1:nGtest){
//         ypred[n] = normal_rng(beta0_t[treetest[n]] + xGtest[n] * u_beta, sigma_y);
//     }
//     //for(n in 1:nGtest){
//         //yrep[n] = gamma_rng(beta0_t[treetest[n]]+xGtest[n]*u_beta,sigma_y);
//         //}
//         
// }


