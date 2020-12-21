
#' lm_linear
#'
#' @param y A vector of length n of tree ring increments on a log scale
#' @param X A n times p matrix of covariate values
#' @param n_plot An integer for the number of plots
#' @param n_tree An integer for the number of trees
#' @param plot_by_tree_idx A vector of length n_tree that indicates which plot each tree was sampled from
#' @param tree_idx A vector of length n indicating which tree the observation came from
#' @param X_pred A n_pred by p matrix of covariate values at which to generate predictions
#' @param tree_idx_pred A vector of length n_pred indicating which tree the observation came from
#' @param ... Additional arguments to `rstan::sampling()`
#'
#' @return A stan object that is the result of fitting the linear regression model
#' @export
#'
lm_linear <- function(y, X, n_plot, n_tree, plot_by_tree_idx, tree_idx, X_pred, tree_idx_pred, ...) {
  ## add in error checking and unit testing

  standata <- list(
    y                = y,
    X                = X,
    n                = length(y),
    K                = ncol(X),
    n_plot           = n_plot,
    n_tree           = n_tree,
    plot_by_tree_idx = plot_by_tree_idx,
    tree_idx         = tree_idx,
    n_pred           = nrow(X_pred),
    X_pred           = X_pred,
    tree_idx_pred    = tree_idx_pred
  )

  out <- rstan::sampling(stanmodels$linear, data = standata, ...)
  return(out)

}


