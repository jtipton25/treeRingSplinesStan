% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_spline_interactions.R
\name{build_spline_interactions}
\alias{build_spline_interactions}
\title{Generate separable B-spline interactions with multivariate regression penalty}
\usage{
build_spline_interactions(X, df)
}
\arguments{
\item{X}{A $n x p$ matrix of covariates from which you want to evaluate
all pairwise combinations of b-spline basis interactions.}

\item{df}{The degree of freedom for the basis function expansion}
}
\value{
An array of separable B-spline interactions. The first dimension is
the number of pairwise interactions (p choose 2). The second dimension is the number
of observations $n$ and the third dimension is the number of knots in the
basis expansion $df^2$.
}
\description{
Generate separable B-spline interactions with multivariate regression penalty
}
\examples{

library(splines)
set.seed(111)
X <- matrix(rnorm(40), 10, 4) ## 10 observations and 4 variables

Z <- build_spline_interactions(X, 5)

}
