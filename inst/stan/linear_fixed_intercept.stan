
data {
    // data
    int<lower=0> K;                        // number of predictor variables
    int<lower=0> n;                        // number of observations
    matrix[n, K] X;                        // matrix of covariates
    vector[n] y;                           // log tree ring increment
    // indexing
    int<lower=0> n_plot;                   // number of plots
    int<lower=0> n_tree;                   // number of trees
    int<lower=1> plot_by_tree_idx[n_tree]; // index for which plot each tree is in
    int<lower=1> tree_idx[n];              // index for trees
    // posterior predictions
    int<lower=0> n_pred;                   // number of posterior predictions
    matrix[n_pred, K] X_pred;              // covariates for posterior predictions
    int<lower=1> tree_idx_pred[n_pred];    // index for trees for posterior predictions
}

parameters {
    // intercept terms
    real mu_beta0;                          // intercept means
    real beta0_p_tilde_fixed[n_plot-1];     // plot-level intercepts (hierarchically centered)
    // real beta0_p_tilde[n_plot];             // plot-level intercepts (hierarchically centered)
    real<lower=0> s_beta0_p;                // plot-level intercept variance (hierarchically centered)
    real beta0_t_tilde_fixed[n_tree-1];     // tree-level intercepts (hierarchically centered)
    // real beta0_t_tilde[n_tree];             // tree-level intercepts (hierarchically centered)
    real<lower=0> s_beta0_t;                // tree-level intercept variance (hierarchically centered)

    // fixed effects
    vector[K] beta;                         // regression fixed effects

    // measurement/data model
    real<lower=0> sigma_y;                  // residual standard deviation

}

transformed parameters {
    real beta0_p[n_plot];                   // plot-level intercepts
    real beta0_t[n_tree];                   // tree-level intercepts
    vector[n] mu;                           // growth model mean response
    real beta0_p_tilde[n_plot];             // plot-level intercepts (hierarchically centered)
    real beta0_t_tilde[n_tree];             // tree-level intercepts (hierarchically centered)


    // fix the last intercept term
    beta0_p_tilde[n_plot] = 0.0;
    for(p in 1:(n_plot - 1)){
        beta0_p_tilde[p] = beta0_p_tilde_fixed[p];
    }
    for(p in 1:n_plot){
        beta0_p[p] = mu_beta0 + s_beta0_p * beta0_p_tilde[p];
    }

    // fix the last intercept term
    beta0_t_tilde[n_tree] = 0.0;
    for(t in 1:(n_tree - 1)){
        beta0_t_tilde[t] = beta0_t_tilde_fixed[t];
    }
    for(t in 1:n_tree){
        beta0_t[t] = beta0_p[plot_by_tree_idx[t]] + s_beta0_t * beta0_t_tilde[t];
    }

    for(i in 1:n) {
        mu[i] = beta0_t[tree_idx[i]] + X[i] * beta;
    }

}

model {

    // priors for intercept terms
    mu_beta0 ~ normal(0, 100);
    beta0_p_tilde_fixed ~ normal(0,1);
    beta0_t_tilde_fixed ~ normal(0,1);
    // beta0_p_tilde ~ normal(0,1);
    // beta0_t_tilde ~ normal(0,1);
    s_beta0_p ~ cauchy(0, 2.5);
    s_beta0_t ~ cauchy(0, 2.5);

    // prior for regression coefficients
    beta ~ normal(0, 100);

    // priors for data/measurement model
    sigma_y ~ cauchy(0, 2.5);

    // likelihood
    y ~ normal(mu, sigma_y);

}

generated quantities{
    vector[n] y_rep;
    vector[n_pred] y_pred;
    vector[n] log_lik;

    for(i in 1:n){
        y_rep[i] = normal_rng(beta0_t[tree_idx[i]] + X[i] * beta, sigma_y);
        log_lik[i] = normal_lpdf(y[i] | beta0_t[tree_idx[i]] + X[i] * beta, sigma_y);
    }
    for(i in 1:n_pred){
        y_pred[i] = normal_rng(beta0_t[tree_idx_pred[i]] + X_pred[i] * beta, sigma_y);
    }

}


