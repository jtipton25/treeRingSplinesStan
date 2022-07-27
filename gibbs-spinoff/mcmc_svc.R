#' Bayesian Spatially Varying Coefficients model
#'
#' this function runs Markov Chain Monte Carlo to estimate the Bayesian spatially-varying coefficient model
#' @param y is a \eqn{n}{n} vector of Gaussian data.
#' @param X is a \eqn{n \times p}{n x p} matrix of fixed effects (like latitude, elevation, etc)
#' @param climate_normal is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param params is a list of parameter settings. The list
#' \code{params} must contain the following values:
#' * \code{n_adapt}: A positive integer number of adaptive MCMC iterations.
#' * \code{n_mcmc}: A positive integer number of total MCMC iterations
#' post adaptation.
#' * \code{n_thin}: A positive integer number of MCMC iterations per saved
#' sample.
#' * \code{n_message}: A positive integer number of frequency of iterations
#'  to output a progress message. For example, \code{n_message = 50}
#'  outputs progress messages every 50 iterations.
#' @param df The degrees of freedom parameter for the basis expansion
#' @param priors is the list of prior settings.
#' @param M The number of resolutions.
#' @param n_neighbors The expected number of neighbors for each interior basis function. This determines the basis radius parameter.
#' @param n_coarse_grid The number of basis functions in one direction (e.g. \code{n_coarse_grid = 10} results in a \eqn{10 \times 10}{10x10} course grid which is further extended by the number of additional padding basis functions given by \code{n_padding}.
#' @param n_padding The number of additional boundary points to add on each boundary. For example, n_padding = 5 will add 5 boundary knots to the both the left  and right side of the grid).
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the initial values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param verbose Should verbose output be printed? Typically this is only useful for troubleshooting.
#' @param joint Should the spatial parameters alpha be sampled jointly or by each resolution
#' @param constraint What constraint should be applied to the spatial process? Options include no constraint (`constraint = "unconstrained"`), a constraint so the entire MRA process sums to 0 (`constraint = "overall"`), a constraint so that each of the M levels of the MRA process sum to 0 (`constraint = "resolution"`), or whether the predicted process must sum to 0 (`constraint = "predicted"`). Note: `constraint = "predicted"` is NOT currently supported.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE} (see spam and Matrix packages for details).
#' @param n_chain is the MCMC chain id. The default is 1.
#'
#'
#' @export
#'
#' @importFrom mvnfast rmvn
#' @importFrom fields fields.rdist.near
#' @importFrom Matrix Cholesky
#' @importFrom stats lm rgamma
#' @importFrom BayesMRA rmvn_arma
#' @import spam
#' @import spam64

mcmc_svc <- function(
  y,
  X,
  climate_normal,
  params,
  df            = 6,
  priors        = NULL,
  inits         = NULL,
  config        = NULL,
  verbose       = FALSE,
  use_spam      = TRUE, ## use spam or Matrix for sparse matrix operations
  n_chain       = 1
) {


  ## add in overall prior term for intercepts -- currently setting these with variance 25


  ##
  ## Run error checks
  ##

  if (!is_numeric_vector(y, length(y)))
    stop("y must be a numeric vector of length N.")
  if (length(y) != nrow(X))
    stop("X must have the same number of rows as the length of y.")
  if (!is_numeric_matrix(X, length(y), ncol(X)))
    stop("X must be a numeric matrix with N rows.")
  if (!is_numeric_matrix(climate_normal, length(y), 2))
    stop("climate_normal must be a numeric matrix with N rows and 2 columns.")

  ## check the params list
  if (!is_positive_integer(params$n_adapt, 1))
    stop("params must contain a positive integer n_adapt.")
  if (!is_positive_integer(params$n_mcmc, 1))
    stop("params must contain a positive integer n_mcmc.")
  if (!is_positive_integer(params$n_thin, 1))
    stop("params must contain a positive integer n_thin.")
  if (!is_positive_integer(params$n_message, 1))
    stop("params must contain a positive integer n_message.")


  params$n_adapt   <- as.integer(params$n_adapt)
  params$n_mcmc    <- as.integer(params$n_mcmc)
  params$n_thin    <- as.integer(params$n_thin)
  params$n_message <- as.integer(params$n_message)

  ## check the priors list


  ## check mcmc input
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("verbose must be either TRUE or FALSE.")
  }

  if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
    stop("use_spam must be either TRUE or FALSE.")
  }

  if (!is_positive_integer(n_chain, 1)) {
    stop("n_chain must be a positive integer")
  }

  n_chain <- as.integer(n_chain)

  ##
  ## setup config
  ##

  ## do we sample the functional relationship parameters? This is primarily
  ## used to troubleshoot model fitting using simulated data

  sample_beta0 <- TRUE
  if (!is.null(config)) {
    if (!is.null(config[['sample_beta0']])) {
      sample_beta0 <- config[['sample_beta0']]
      if (!is.logical(sample_beta0) | is.na(sample_beta0))
        stop('If specified, sample_beta0 must be TRUE or FALSE')
    }
  }

  sample_beta <- TRUE
  if (!is.null(config)) {
    if (!is.null(config[['sample_beta']])) {
      sample_beta <- config[['sample_beta']]
      if (!is.logical(sample_beta) | is.na(sample_beta))
        stop('If specified, sample_beta must be TRUE or FALSE')
    }
  }


  # sample_rho <- TRUE
  # if (!is.null(config)) {
  #     if (!is.null(config[['sample_rho']])) {
  #         sample_rho <- config[['sample_rho']]
  #     }
  # }

  ## do we sample the regression variance parameters? This is primarily
  ## used to troubleshoot model fitting using simulated data
  sample_tau20 <- TRUE
  if (!is.null(config)) {
    if (!is.null(config[['sample_tau20']])) {
      sample_tau20 <- config[['sample_tau20']]
      if (!is.logical(sample_tau20) | is.na(sample_tau20))
        stop('If specified, sample_tau20 must be TRUE or FALSE')
    }
  }

  ## do we sample the regression variance parameters? This is primarily
  ## used to troubleshoot model fitting using simulated data
  sample_tau2 <- TRUE
  if (!is.null(config)) {
    if (!is.null(config[['sample_tau2']])) {
      sample_tau2 <- config[['sample_tau2']]
      if (!is.logical(sample_tau2) | is.na(sample_tau2))
        stop('If specified, sample_tau2 must be TRUE or FALSE')
    }
  }

  ## do we sample the nugger variance parameter
  sample_sigma2 <- TRUE
  if (!is.null(config)) {
    if (!is.null(config[['sample_sigma2']])) {
      sample_sigma2 <- config[['sample_sigma2']]
      if (!is.logical(sample_sigma2) | is.na(sample_sigma2))
        stop('If specified, sample_sigma2 must be TRUE or FALSE')
    }
  }

  N      <- length(y)
  # n_time <- dim(Y)[3]
  p      <- ncol(X)

  ##
  ## center and scale the input and covariates
  ##

  sd_y <- sd(y)
  mu_y <- mean(y)
  # y <- (y - mu_y) / sd_y

  mu_X <- apply(X[, -1, drop = FALSE], 2, mean)
  sd_X <- apply(X[, -1, drop = FALSE], 2, sd)
  # for(i in 2:ncol(X)) {
  #   X[, i] <- (X[, i] - mu_X[i-1]) / sd_X[i-1]
  # }


  ##
  ## setup climate normal basis
  ##

  W <- cbind(1, build_spline_interactions(cbind(climate_normal[, 1], climate_normal[, 2]), df = df)[1, , ])
  W <- as.spam(W)

  ##
  ## initial values
  ##

  tX  <- t(X)
  tXX <- tX %*% X
  tW <- t(W)
  tWW <- tW %*% W
  tXW   <- list(length = p)
  tXWWX <- list(length = p)
  for (j in 1:p) {
    tXW[[j]]   <- as.spam(t(X[, j] * as.matrix(W)))
    tXWWX[[j]] <- tXW[[j]] %*% as.spam(X[, j] * as.matrix(W))
  }

  ##
  ## priors for beta
  ##

  mu_beta        <- rep(0, p)
  Sigma_beta     <- 10 * diag(p)

  # ## check if priors for mu_beta are specified
  # if (!is.null(priors[['mu_beta']])) {
  #   if(!is_numeric_vector(priors[['mu_beta']], p))
  #     stop("If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
  #   if (all(!is.na(priors[['mu_beta']]))) {
  #     mu_beta <- priors[['mu_beta']]
  #   }
  # }
  #
  # ## check if priors for Sigma_beta are specified
  # if (!is.null(priors[['Sigma_beta']])) {
  #   if(!is_sympd_matrix(priors[['Sigma_beta']], p))
  #     stop("If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
  #   if (all(!is.na(priors[['Sigma_beta']]))) {
  #     Sigma_beta <- priors[['Sigma_beta']]
  #   }
  # }
  # Sigma_beta_chol <- tryCatch(
  #   chol(Sigma_beta),
  #   error = function(e) {
  #     if (verbose)
  #       message("The Cholesky decomposition of the prior covariance Sigma_beta was ill-conditioned and mildy regularized.")
  #     chol(Sigma_beta + 1e-8 * diag(p))
  #   }
  # )
  # Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)

  ##
  ## initialize sigma2
  ##

  alpha_sigma2 <- 1
  beta_sigma2  <- 1
  ## check if priors for alpha_sigma2 are specified
  if (!is.null(priors[['alpha_sigma2']])) {
    if (!is_positive_numeric(priors[['alpha_sigma2']], 1))
      stop("If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    if (all(!is.na(priors[['alpha_sigma2']]))) {
      alpha_sigma2 <- priors[['alpha_sigma2']]
    }
  }
  ## check if priors for beta_sigma2 are specified
  if (!is.null(priors[['beta_sigma2']])) {
    if (!is_positive_numeric(priors[['beta_sigma2']], 1))
      stop("If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    if (all(!is.na(priors[['beta_sigma2']]))) {
      beta_sigma2 <- priors[['beta_sigma2']]
    }
  }

  sigma2  <- pmax(pmin(1 / rgamma(1, alpha_sigma2, beta_sigma2), 5), 0.1)
  sigma   <- sqrt(sigma2)

  ##
  ## Initialize tau2
  ##
  tau20 <- 1
  tau2   <- rep(1, p)

  alpha_tau2 <- 0.01
  beta_tau2  <- 0.01
  ## check if priors for alpha_tau2 are specified
  if (!is.null(priors[['alpha_tau2']])) {
    if (!is_positive_numeric(priors[['alpha_tau2']], 1))
      stop("If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    if (all(!is.na(priors[['alpha_tau2']]))) {
      alpha_tau2 <- priors[['alpha_tau2']]
    }
  }
  ## check if priors for beta_tau2 are specified
  if (!is.null(priors[['beta_tau2']])) {
    if (!is_positive_numeric(priors[['beta_tau2']], 1))
      stop("If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    if (all(!is.na(priors[['beta_tau2']]))) {
      beta_tau2 <- priors[['beta_tau2']]
    }
  }

  ##
  ## initialize betas
  ##

  Q0 <- bdiag.spam(1 / 25, tau20 * make_Q(df, 0.99))
  Q <- list(length = p)
  for (j in 1:p) {
    Q[[j]] <- bdiag.spam(1 / 25, tau2[j] * make_Q(df, 0.99))
  }
  beta0 <- drop(rmvnorm.prec(1, rep(0, df^2+1), Q0))
  Rstruct0 <- chol(1 / sigma2 * tWW + Q0)

  beta <- matrix(0, df^2 + 1, p)
  for (j in 1:p) {
    beta[, j] <- rmvnorm.prec(1, rep(0, df^2+1), Q[[j]])
  }
  Rstruct <- list(length = p)
  for (j in 1:p) {
    Rstruct[[j]] <- chol(1 / sigma2 * tXWWX[[j]] + Q[[j]])
  }

  Wbeta0 <- W %*% beta0
  Wbeta  <- W %*% beta
  XWbeta <- X * Wbeta

  ## setup sparse cholesky structures -- is this worth it?
  ## check timings here to see if this is worthwhile
#
#   A0       <- 1 / sigma2 * tWW + Q0
#   b0       <- 1 / sigma2 * tW %*% (y - apply(XWbeta, 1, sum))
#   Rstruct0 <- chol(A0)



  ##
  ## intialize an ICAR structure to initialize the parameter alpha
  ##
  Q0 <- bdiag.spam(1 / 25, tau20 * make_Q(df, 1))
  Q <- list(length = p)
  for (j in 1:p) {
    Q[[j]] <- bdiag.spam(1 / 25, tau2[j] * make_Q(df, 1))
  }

  ##
  ## sampler config options -- to be added later
  ##
  #

  ## check for initial values
  ##

  ## initial values for beta0
  if (!is.null(inits[['beta0']])) {
    if(!is_numeric_vector(inits[['beta0']], p))
      stop("initial value for beta0 must be a numeric vector of length p")
    if (all(!is.na(inits[['beta0']]))) {
      beta0 <- inits[['beta0']]
    }
  }
  Wbeta0 <- W %*% beta0

  ## initial values for beta
  if (!is.null(inits[['beta']])) {
    if(!is_numeric_vector(inits[['beta']], p))
      stop("initial value for beta must be a numeric vector of length p")
    if (all(!is.na(inits[['beta']]))) {
      beta <- inits[['beta']]
    }
  }
  Wbeta  <- W %*% beta
  XWbeta <- X * Wbeta

  ## intial values for sigma2
  if (!is.null(inits[['sigma2']])) {
    if(!is_positive_numeric(inits[['sigma2']], 1))
      stop("initial value for sigma2 must be a positive numeric value")
    if (all(!is.na(inits[['sigma2']]))) {
      sigma2 <- inits[['sigma2']]
    }
  }

  ## initial values for tau20
  if (!is.null(inits[['tau20']])) {
    if(!is_positive_numeric(inits[['tau20']], M) | !is.vector(inits[['tau20']]))
      stop("initial value for tau20 must be a positive numeric vector of length M")
    if (all(!is.na(inits[['tau20']]))) {
      tau20 <- inits[['tau20']]
    }
  }
  # update Q0 with the initial value
  Q0 <- bdiag.spam(1 / 25, tau20 * make_Q(df, 1))

  ## initial values for tau2
  if (!is.null(inits[['tau2']])) {
    if(!is_positive_numeric(inits[['tau2']], M) | !is.vector(inits[['tau2']]))
      stop("initial value for tau2 must be a positive numeric vector of length M")
    if (all(!is.na(inits[['tau2']]))) {
      tau2 <- inits[['tau2']]
    }
  }
  # update Q with the initial value
  for (j in 1:p) {
    Q[[j]] <- bdiag.spam(1 / 25, tau2[j] * make_Q(df, 1))
  }


  ##
  ## setup save variables
  ##

  n_save       <- params$n_mcmc / params$n_thin
  beta0_save   <- matrix(0, n_save, df^2 + 1)
  beta_save    <- array(0, dim = c(n_save, df^2 + 1, p))
  tau20_save   <- rep(0, n_save)
  tau2_save    <- matrix(0, n_save, p)
  sigma2_save  <- rep(0, n_save)

  ##
  ## initialize tuning
  ##

  ##
  ## tuning variables for adaptive MCMC
  ##


  ##
  ## Starting MCMC chain
  ##

  message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations \n")

  for (k in 1:(params$n_adapt + params$n_mcmc)) {
    if (k == params$n_adapt + 1) {
      message("Starting MCMC fitting for chain ", n_chain, ", running for ", params$n_mcmc, " iterations \n")
    }
    if (k %% params$n_message == 0) {
      if (k <= params$n_adapt) {
        message("MCMC adaptation iteration ", k, " for chain ", n_chain)
      } else {
        message("MCMC fitting iteration ", k - params$n_adapt, " for chain ", n_chain)
      }
    }

    ##
    ## sample sigma2
    ##

    if (sample_sigma2) {
      if (verbose)
        message("sample sigma2")

      devs   <- y - Wbeta0 - apply(XWbeta, 1, sum)
      SS     <- sum(devs^2)
      sigma2 <- 1 / rgamma(1, N / 2 + alpha_sigma2, SS / 2 + beta_sigma2)
      sigma  <- sqrt(sigma2)
    }

    ##
    ## sample beta0 -- double check these values
    ##

    if (sample_beta0) {
      if (verbose)
        message("sample beta0")

      # A      <- as.matrix(1 / sigma2 * tWW + Q0)
      A     <- 1 / sigma2 * tWW + Q0
      b     <- 1 / sigma2 * tW %*% (y - apply(XWbeta, 1, sum))
      ## guarantee a symmetric matrix
      A     <- (A + t(A)) / 2
      # beta0  <- drop(rmvn_arma(A, b))
      beta0 <- drop(rmvnorm.canonical(1, b, A, Rstruct = Rstruct0))

      ## update Wbeta0
      Wbeta0 <- drop(W %*% beta0)
    }

    ##
    ## sample beta -- double check these values
    ##

    if (sample_beta) {
      if (verbose)
        message("sample beta")
      # A <- 1 / sigma2 * tXWXW[[j]] + Q
      # b <- 1 / sigma1 * tW(y - W %*% beta0 - apply(X[, -j] * (W %*% beta)[, -j], 1, sum))

      for (j in 1:p) {

        # A <- 1 / sigma2 * tXWWX[, , j] + Q[[j]]
        A <- 1 / sigma2 * tXWWX[[j]] + Q[[j]]
        # b <- drop(1 / sigma2 * tXW[, , j] %*% (y - Wbeta0 - apply(XWbeta[, -j, drop = FALSE], 1, sum)))
        b <- drop(1 / sigma2 * tXW[[j]] %*% (y - Wbeta0 - apply(XWbeta[, -j, drop = FALSE], 1, sum)))
        ## guarantee a symmetric matrix
        A         <- (A + t(A)) / 2
        # beta[, j] <- drop(rmvn_arma(A, b))
        beta[, j] <- drop(rmvnorm.canonical(1, b, A, Rstruct = Rstruct[[j]]))

        ## update XWbeta
        Wbeta[, j] <- W %*% beta[, j]
        XWbeta[, j] <- X[, j] * Wbeta[, j]
      }
    }

    ##
    ## sample tau20
    ##

    if (sample_tau20) {
      if (verbose)
        message("sample tau20")
      devs  <- beta0[-1]
      SS    <- as.numeric(devs %*% (Q0[-1, -1] %*% devs))
      tau20 <- rgamma(1, alpha_tau2 + df^2 / 2, beta_tau2 + SS / 2)
    }

    Q0 <- bdiag.spam(1 / 25, tau20 * make_Q(df, 1))


    ##
    ## sample tau2
    ##

    if (sample_tau2) {
      if (verbose)
        message("sample tau2")
      for (j in 1:p) {
        devs    <- beta[-1, j]
        SS      <- as.numeric(devs %*% (Q[[j]][-1, -1] %*% devs))
        tau2[j] <- rgamma(1, alpha_tau2 + df^2 / 2, beta_tau2 + SS / 2)
        Q[[j]] <- bdiag.spam(1 / 25, tau2[j] * make_Q(df, 1))
      }
    }


    ##
    ## save variables
    ##

    if (k >= params$n_adapt) {
      if (k %% params$n_thin == 0) {
        save_idx                <- (k - params$n_adapt) / params$n_thin
        beta0_save[save_idx, ]  <- beta0
        beta_save[save_idx, , ] <- beta
        tau20_save[save_idx]    <- tau20
        tau2_save[save_idx, ]   <- tau2
        sigma2_save[save_idx]   <- sigma2
      }
    }

    ##
    ## End of MCMC loop
    ##
  }

  ## print out acceptance rates -- no tuning in this model


  ##
  ## return the MCMC output -- think about a better way to make this a class
  ##

  # for(i in 2:ncol(X)) {
  #   X[, i] <- X[, i] * sd_X[i-1] + mu_X[i-1]
  # }

  out <- list(
    beta0  = beta0_save,
    beta   = beta_save,
    tau20  = tau20_save,
    tau2   = tau2_save,
    sigma2 = sigma2_save,
    W      = W,
    data   = list(
      y    = y * sd_y + mu_y,
      mu_y = mu_y,
      sd_y = sd_y,
      X    = X,
      mu_X = mu_X,
      sd_X = sd_X,
      climate_normal = climate_normal),
    model  = list(
      params     = params,
      priors     = priors,
      inits      = inits,
      config     = config,
      verbose    = verbose,
      use_spam   = use_spam,
      n_chain    = n_chain)
    ## add in run-time variables as an object too
  )

  class(out) <- "mcmc_svc"

  return(out)
}
