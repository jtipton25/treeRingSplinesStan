
data {
    // data
    int<lower=0> K;                        // number of predictor variables
    int<lower=0> n;                        // number of observations
    int<lower=1> n_tree;                   // number of trees
    matrix[n, n_tree] X_int;               // matrix of intercept covariates
    matrix[n, K] X;                        // matrix of covariates
    vector[n] y;                           // log tree ring increment
    // indexing
    // posterior predictions
    int<lower=0> n_pred;                   // number of posterior predictions
    matrix[n_pred, K] X_pred;              // covariates for posterior predictions
    matrix[n_pred, n_tree] X_int_pred;     // covariates for posterior predictions
}

parameters {
    // intercept terms
    vector[n_tree] beta_int;                // intercept coefficients

    // fixed effects
    vector[K] beta;                         // regression fixed effects

    // measurement/data model
    real<lower=0> sigma_y;                  // residual standard deviation

}


model {

    // priors for intercept terms
    beta_int ~ normal(0, 10);

    // prior for regression coefficients
    beta ~ normal(0, 10);

    // priors for data/measurement model
    sigma_y ~ cauchy(0, 2.5);

    // likelihood
    y ~ normal(X_int * beta_int + X * beta, sigma_y);

}

generated quantities{
    vector[n] y_rep;
    vector[n_pred] y_pred;
    vector[n] log_lik;

    for(i in 1:n){
        y_rep[i] = normal_rng(X_int[i] * beta_int + X[i] * beta, sigma_y);
        log_lik[i] = normal_lpdf(y[i] | X_int[i] * beta_int + X[i] * beta, sigma_y);
    }
    for(i in 1:n_pred){
        y_pred[i] = normal_rng(X_int_pred[i] * beta_int + X_pred[i] * beta, sigma_y);
    }

}


