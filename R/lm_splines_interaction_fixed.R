
#' lm_splines_interaction_fixed
#'
#' @param y A vector of length n of tree ring increments on a log scale
#' @param X A n times p matrix of covariate values
#' @param X_bs A p times n times df array of covariate values bspline basis function expansions
#' @param X_bs_int A p times n times df_int array of covariate values bspline basis function expansions
#' @param Q A df_int times df_int matrix of spline penalty parameters
#' @param n_plot An integer for the number of plots
#' @param n_tree An integer for the number of trees
#' @param plot_by_tree_idx A vector of length n_tree that indicates which plot each tree was sampled from
#' @param tree_idx A vector of length n indicating which tree the observation came from
#' @param X_pred A n_pred by p matrix of covariate values at which to generate predictions
#' @param X_bs_pred A p times n_pred times df matrix of covariate values at which to generate predictions
#' @param X_bs_int_pred A n_interaction times n_pred times df^2 matrix of covariate values at which to generate predictions
#' @param tree_idx_pred A vector of length n_pred indicating which tree the observation came from
#' @param ... Additional argument to `rstan::sampling()`
#'
#' @return A stan object that is the result of fitting the spline interaction regression model
#' @export
#'
lm_splines_interaction_fixed <- function(y, X, X_bs, X_bs_int, Q, n_plot, n_tree, plot_by_tree_idx, tree_idx, X_pred, X_bs_pred, X_bs_int_pred, tree_idx_pred, ...) {
  ## add in error checking and unit testing

if (nrow(Q) != dim(X_bs_int)[3])
  stop("The b-spline interaction matrix must have the same third dimension size as the penalty matrix Q")

  standata <- list(
    y                = y,
    X                = X,
    X_bs             = X_bs,
    X_bs_int         = X_bs_int,
    Q                = Q,
    q                = nrow(X_bs_int),
    df_int           = nrow(Q),
    p                = dim(X_bs)[1],
    df               = dim(X_bs)[3],
    n                = length(y),
    K                = ncol(X),
    n_plot           = n_plot,
    n_tree           = n_tree,
    plot_by_tree_idx = plot_by_tree_idx,
    tree_idx         = tree_idx,
    n_pred           = nrow(X_pred),
    X_pred           = X_pred,
    X_bs_pred        = X_bs_pred,
    X_bs_int_pred    = X_bs_int_pred,
    tree_idx_pred.   = tree_idx_pred
  )

  out <- rstan::sampling(stanmodels$splines_interaction_fixed, data = standata, ...)
  return(out)

}


