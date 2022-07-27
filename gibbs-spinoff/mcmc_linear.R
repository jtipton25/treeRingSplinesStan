
#' lm_linear
#'
#' @param y A vector of length n of tree ring increments on a log scale
#' @param X A n times p matrix of covariate values
#' @param plot_idx A vector of length n that indicates which plot each observation was sampled from
#' @param tree_idx A vector of length n indicating which tree the observation came from
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
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param verbose Should verbose output be printed? Typically this is only useful for troubleshooting.
#' @param n_chain is the MCMC chain id. The default is 1.
#'
#' @return A MCMC object that is the result of fitting the linear regression model
#' @importFrom BayesMRA rmvn_arma
#' @export
#'
mcmc_linear <- function(y, X, plot_idx, tree_idx, params, config = NULL, verbose = FALSE, n_chain = 1) {

  ## check the params list
  # if (!is_positive_integer(params$n_adapt, 1))
  #   stop("params must contain a positive integer n_adapt.")
  # if (!is_positive_integer(params$n_mcmc, 1))
  #   stop("params must contain a positive integer n_mcmc.")
  # if (!is_positive_integer(params$n_thin, 1))
  #   stop("params must contain a positive integer n_thin.")
  # if (!is_positive_integer(params$n_message, 1))
  #   stop("params must contain a positive integer n_message.")

  ## do we sample the residual error variance?
  sample_sigma2 <- TRUE
  if (!is.null(config)) {
    if (!is.null(config[['sample_sigma2']])) {
      sample_sigma2 <- config[['sample_sigma2']]
      if (!is.logical(sample_sigma2) | is.na(sample_sigma2))
        stop('If specified, sample_sigma2 must be TRUE or FALSE')
    }
  }

  ## do we sample the intercept parameters?
  sample_beta_intercept <- TRUE
  if (!is.null(config)) {
    if (!is.null(config[['sample_beta_intercept']])) {
      sample_beta_intercept <- config[['sample_beta_intercept']]
      if (!is.logical(sample_beta_intercept) | is.na(sample_beta_intercept))
        stop('If specified, sample_beta_intercept must be TRUE or FALSE')
    }
  }
  ## do we sample the slope parameters?
  sample_beta <- TRUE
  if (!is.null(config)) {
    if (!is.null(config[['sample_beta']])) {
      sample_beta <- config[['sample_beta']]
      if (!is.logical(sample_beta) | is.na(sample_beta))
        stop('If specified, sample_beta must be TRUE or FALSE')
    }
  }

  n_adapt   <- as.integer(params$n_adapt)
  n_mcmc    <- as.integer(params$n_mcmc)
  n_thin    <- as.integer(params$n_thin)
  n_message <- as.integer(params$n_message)


  ## add in error checking and unit testing



  # precalculate
  N <- length(y)
  tX <- t(X)
  tXX <- tX %*% X

  X_intercept <- model.matrix(~ factor(plot_idx) + factor(tree_idx))
  intercept_idx <- c("intercept", rep("plot", n_plot - 1), rep("tree", n_tree - 1))
  n_intercept <- c(
    sum(intercept_idx == "intercept"),
    sum(intercept_idx == "plot"),
    sum(intercept_idx == "tree"))
  tX_int <- t(X_intercept)
  tXX_int <- tX_int %*% X_intercept

  mu_beta_intercept <- rep(0, ncol(X_intercept))
  ## vague prior on the scale of the data
  s2_beta_p <- runif(1, 0.1, 1)
  s2_beta_t <- runif(1, 0.1, 1)
  Sigma_beta_intercept <- diag(5, rep(s2_beta_p, n_intercept[2]), rep(s2_beta_t, n_intercept[3]))
  Q_beta_intercept <-  diag(1 / c(5, rep(s2_beta_p, n_intercept[2]), rep(s2_beta_t, n_intercept[3])))

  # initialize the intercept term
  beta_intercept <- rnorm(ncol(X_intercept))
  # initialize the slope terms
  beta <- rnorm(ncol(X))
  mu_beta <- rep(0, ncol(X))
  ## vague prior on the scale of the data
  Sigma_beta <- 5 * diag(ncol(X))
  Q_beta <- 1 / 5 * diag(ncol(X))

  ## intialize the residual error
  alpha_sigma2 <- 0.5
  beta_sigma2 <- 0.5

  sigma2 <- runif(1, 0.1, 1)

  n_save              <- n_mcmc / n_thin
  beta_intercept_save <- matrix(0, n_save, ncol(X_intercept))
  s2_beta_p_save      <- rep(0, n_save)
  s2_beta_t_save      <- rep(0, n_save)
  beta_save           <- matrix(0, n_save, ncol(X))
  sigma2_save         <- rep(0, n_save)
  lp_save             <- matrix(0, n_save, N)

  message("Starting MCMC for chain ", n_chain, ", running for ", n_adapt, " adaptive iterations and ", n_mcmc, " fitting iterations \n")

  for (k in 1:(n_adapt + n_mcmc)) {

    if (k == n_adapt + 1) {
      message("Starting MCMC fitting for chain ", n_chain, ", running for ", n_mcmc, " iterations \n")
    }
    if (k %% n_message == 0) {
      if (k <= n_adapt) {
        message("MCMC adaptation iteration ", k, " for chain ", n_chain)
      } else {
        message("MCMC fitting iteration ", k - n_adapt, " for chain ", n_chain)
      }
    }

    ## sample beta_intercept
    if (sample_beta_intercept) {
      if (verbose)
        message("sample beta_intercept")

      A              <- 1 / sigma2 * tXX_int + Q_beta_intercept
      b              <- 1 / sigma2 * tX_int %*% (y - X %*% beta) + Q_beta_intercept %*% mu_beta_intercept
      A              <- (A + t(A)) / 2 # guarantee a symmetric matrix
      beta_intercept <- rmvn_arma(A, b)
    }

    ## sample beta slopes
    if (sample_beta) {
      if (verbose)
        message("sample beta")

      A      <- 1 / sigma2 * tXX + Q_beta
      b      <- 1 / sigma2 * tX %*% (y - X_intercept %*% beta_intercept) + Q_beta %*% mu_beta
      A      <- (A + t(A)) / 2 # guarantee a symmetric matrix
      beta   <- rmvn_arma(A, b)
    }

    ## sample sigma2

    if (sample_sigma2) {
      if (verbose)
        message("sample sigma2")

      devs <- y - X_intercept %*% beta_intercept - X %*% beta
      SS        <- sum(devs^2)
      sigma2 <- 1 / rgamma(1, N / 2 + alpha_sigma2, SS / 2 + beta_sigma2)
      sigma     <- sqrt(sigma2)
    }

    ## sample s_beta_p

    devs      <- beta_intercept[intercept_idx == "plot"]
    SS        <- sum(devs)^2
    s2_beta_p <- 1 / rgamma(1, n_intercept[2] / 2 + 0.5, SS / 2 + 0.5)

    ## sample_s_beta_t
    devs      <- beta_intercept[intercept_idx == "tree"]
    SS        <- sum(devs)^2
    s2_beta_t <- 1 / rgamma(1, n_intercept[3] / 2 + 0.5, SS / 2 + 0.5)

    Q_beta_intercept <-  diag(1 / c(5, rep(s2_beta_p, n_intercept[2]), rep(s2_beta_t, n_intercept[3])))

    ##
    ## save variables
    ##

    if (k >= n_adapt) {
      if (k %% n_thin == 0) {
        save_idx                <- (k - n_adapt) / n_thin
        beta_intercept_save[save_idx, ] <- beta_intercept
        beta_save[save_idx, ]           <- beta
        sigma2_save[save_idx]           <- sigma2
        s2_beta_p_save[save_idx]        <- s2_beta_p
        s2_beta_t_save[save_idx]        <- s2_beta_t
        lp_save[save_idx, ]             <- dnorm(y, X_intercept %*% beta_intercept + X %*% beta, sqrt(sigma2), log = TRUE)
      }
    }
  }

  ##
  ## return the MCMC output -- think about a better way to make this a class
  ##

  out <- list(
    beta           = beta_save,
    beta_intercept = beta_intercept_save,
    sigma2         = sigma2_save,
    s2_beta_p      = s2_beta_p_save,
    s2_beta_t      = s2_beta_t_save,
    lp__           = lp_save
  )

  class(out) <- "mcmc_linear"

  return(out)
}
