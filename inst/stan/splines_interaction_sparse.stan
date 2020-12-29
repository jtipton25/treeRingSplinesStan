
functions {
  /**
  * Return the log probability of a proper intrinsic autoregressive (IAR) prior
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a IAR prior
  * @param tau Precision parameter for the IAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of IAR prior up to additive constant
  */
  real sparse_iar_lpdf(vector phi, real tau,
    int[,] W_sparse, vector D_sparse, vector lambda, int df_int, int W_n) {
      row_vector[df_int] phit_D; // phi' * D
      row_vector[df_int] phit_W; // phi' * W
      vector[df_int] ldet_terms;

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, df_int);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      return 0.5 * ((df_int-1) * log(tau)
                    - tau * (phit_D * phi - (phit_W * phi)));
  }
}

data {
  // data
  int<lower=0> n;                          // number of observations
  int<lower=0> K;                          // number of linear predictor variables
  int<lower=0> p;                          // number of non-linear (spline) predictor variables
  int<lower=0> df;                         // number of spline basis functions per predictor variable
  matrix[n, K] X;                          // matrix of covariates with linear responses
  matrix[n, df] X_bs[p];                   // array of matrices of covariates with non-linear (spline) responses
  int<lower=0> q;                          // number of non-linear (spline) predictor variables
  int<lower=0> df_int;                     // number of spline basis functions per predictor variable
  matrix[n, df_int] X_bs_int[q];           // array of matrices of covariates with non-linear (spline) responses
  matrix[df_int, df_int] W;                // spline interaction penalty adjacency matrix
  int W_n;                                 // number of adjacent pairs
  vector[n] y;                             // log tree ring increment
  // indexing
  int<lower=0> n_plot;                     // number of plots
  int<lower=0> n_tree;                     // number of trees
  int<lower=1> plot_by_tree_idx[n_tree];   // index for which plot each tree is in
  int<lower=1> tree_idx[n];                // index for trees
  // posterior predictions
  int<lower=0> n_pred;                     // number of posterior predictions
  matrix[n_pred, K] X_pred;                // covariates for posterior predictions
  matrix[n_pred, df] X_bs_pred[p];         // array of matrices of covariates with non-linear (spline) responses for posterior predictions
  matrix[n_pred, df_int] X_bs_int_pred[q]; // array of matrices of covariates with non-linear (spline) responses for posterior predictions
  int<lower=1> tree_idx_pred[n_pred];      // index for trees for posterior predictions
}

transformed data {
  int W_sparse[W_n, 2];                   // adjacency pairs
  vector[df_int] D_sparse;                // diagonal of D (number of neigbors for each site)
  vector[df_int] lambda;                 // eigenvalues of invsqrtD * W * invsqrtD
  vector[df_int] zeros;

  for (j in 1:df_int) {
    zeros[j] = 0.0;
  }

  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(df_int - 1)) {
      for (j in (i + 1):df_int) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:df_int) D_sparse[i] = sum(W[i]);
  {
    vector[df_int] invsqrtD;
    for (i in 1:df_int) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }

}
parameters {
  // intercept terms
  real mu_beta0;                          // intercept means
  real beta0_p_tilde[n_plot];             // plot-level intercepts (hierarchically centered)
  real<lower=0> s_beta0_p;                // plot-level intercept variance (hierarchically centered)
  real beta0_t_tilde[n_tree];             // tree-level intercepts (hierarchically centered)
  real<lower=0> s_beta0_t;                // tree-level intercept variance (hierarchically centered)

  // fixed effects
  vector[K] beta;                         // regression fixed effects

  // spline (non-linear) effects
  vector[df] beta_bs[p];                  // regression mixed (non-linear/spline) effects
  real<lower=0> sigma_beta[p];            // spline coefficient variances

  // interaction spline (non-linear) effects
  vector[df_int] beta_int[q];             // regression mixed (non-linear/spline) effects
  real<lower=0> sigma_beta_int[q];        // spline coefficient variances

  // measurement/data model
  real<lower=0> sigma_y;                  // residual standard deviation

}

transformed parameters {
  real beta0_p[n_plot];                   // plot-level intercepts
  real beta0_t[n_tree];                   // tree-level intercepts
  vector[n] mu;                           // growth model mean response
  real<lower=0> tau_beta_int[q];          // spline coefficient variances

  for(j in 1:n_plot){
    beta0_p[j] = mu_beta0 + s_beta0_p * beta0_p_tilde[j];
  }

  for(t in 1:n_tree){
    beta0_t[t] = beta0_p[plot_by_tree_idx[t]] + s_beta0_t * beta0_t_tilde[t];
  }

  for(i in 1:n) {
    // temporary variable for mixed effects (non-linear/spline)
    real tmp = 0;
    real tmp_int = 0;

    // spline effects
    for (j in 1:p) {
      tmp = tmp + X_bs[j][i] * beta_bs[j];
    }
    // spline interaction effects
    for (j in 1:q) {
      tmp_int = tmp_int + X_bs_int[j][i] * beta_int[j];
    }
    mu[i] = beta0_t[tree_idx[i]] + X[i] * beta + tmp + tmp_int;
  }

  for (j in 1:q) {
    tau_beta_int[j] = 1 / sigma_beta_int[j];
  }
}

model {

  // priors for intercept terms
  mu_beta0 ~ normal(0, 100);
  beta0_p_tilde ~ normal(0,1);
  beta0_t_tilde ~ normal(0,1);
  s_beta0_p ~ cauchy(0, 2.5);
  s_beta0_t ~ cauchy(0, 2.5);

  // prior for regression coefficients
  beta ~ normal(0, 100);

  // prior for mixed effects (non-linear/spline) coefficients
  sigma_beta ~ cauchy(0, 2.5);
  for (i in 1:p) {
    beta_bs[i][1] ~ normal(0, 1);
    for (k in 2:df) {
      beta_bs[i][k] ~ normal(beta_bs[i][k-1], sigma_beta[i]);
    }
  }

  // prior for interaction effects
  sigma_beta_int ~ cauchy(0, 2.5);
  for (j in 1:q) {
    beta_int[j] ~ sparse_iar(tau_beta_int[j], W_sparse, D_sparse, lambda, df_int, W_n);
  }

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
    // temporary variable for mixed effects (non-linear/spline)
    real tmp = 0;
    real tmp_int = 0;
    for (j in 1:p) {
      tmp = tmp + X_bs[j][i] * beta_bs[j];
    }
    for (j in 1:q) {
      tmp_int = tmp_int + X_bs_int[j][i] * beta_int[j];
    }

    y_rep[i] = normal_rng(beta0_t[tree_idx[i]] + X[i] * beta + tmp, sigma_y);
    log_lik[i] = normal_lpdf(y[i] | beta0_t[tree_idx[i]] + X[i] * beta + tmp, sigma_y);
  }
  for(i in 1:n_pred){
    // temporary variable for mixed effects (non-linear/spline)
    real tmp = 0;
    real tmp_int = 0;
    for (j in 1:p) {
      tmp = tmp + X_bs_pred[j][i] * beta_bs[j];
    }
    for (j in 1:q) {
      tmp_int = tmp_int + X_bs_int_pred[j][i] * beta_int[j];
    }
    y_pred[i] = normal_rng(beta0_t[tree_idx_pred[i]] + X_pred[i] * beta + tmp, sigma_y);
  }

}


