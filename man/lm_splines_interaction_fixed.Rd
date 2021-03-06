% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lm_splines_interaction_fixed.R
\name{lm_splines_interaction_fixed}
\alias{lm_splines_interaction_fixed}
\title{lm_splines_interaction_fixed}
\usage{
lm_splines_interaction_fixed(
  y,
  X,
  X_bs,
  X_bs_int,
  Q,
  n_plot,
  n_tree,
  plot_by_tree_idx,
  tree_idx,
  X_pred,
  X_bs_pred,
  X_bs_int_pred,
  tree_idx_pred,
  ...
)
}
\arguments{
\item{y}{A vector of length n of tree ring increments on a log scale}

\item{X}{A n times p matrix of covariate values}

\item{X_bs}{A p times n times df array of covariate values bspline basis function expansions}

\item{X_bs_int}{A p times n times df_int array of covariate values bspline basis function expansions}

\item{Q}{A df_int times df_int matrix of spline penalty parameters}

\item{n_plot}{An integer for the number of plots}

\item{n_tree}{An integer for the number of trees}

\item{plot_by_tree_idx}{A vector of length n_tree that indicates which plot each tree was sampled from}

\item{tree_idx}{A vector of length n indicating which tree the observation came from}

\item{X_pred}{A n_pred by p matrix of covariate values at which to generate predictions}

\item{X_bs_pred}{A p times n_pred times df matrix of covariate values at which to generate predictions}

\item{X_bs_int_pred}{A n_interaction times n_pred times df^2 matrix of covariate values at which to generate predictions}

\item{tree_idx_pred}{A vector of length n_pred indicating which tree the observation came from}

\item{...}{Additional argument to \code{rstan::sampling()}}
}
\value{
A stan object that is the result of fitting the spline interaction regression model
}
\description{
lm_splines_interaction_fixed
}
