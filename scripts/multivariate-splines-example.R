library(treeRingSplinesStan)
library(igraph)
library(ggplot2)
library(tidyverse)
library(mvnfast)

set.seed(11)
n <- 20
x1 <- seq(0, 1, length = n)
x2 <- seq(0, 1, length = n)
covariates <- expand.grid(x1, x2)

df      <- 12
Q       <- spline_interactions_penalty(df)
idx     <- spline_interactions_lattice(df)

# check the indxing vs. number of neighbors
# cbind(idx, n = apply(as_adjacency_matrix(make_lattice(length = df, dim = 2)), 2, sum))

X_bs    <- spline_interactions(covariates[, 1], covariates[, 2], df)

dat <- data.frame(
  spline      = c(X_bs),
  observation = factor(rep(1:n^2, times = df^2)),
  knot_x      = factor(rep(idx[, 1], each = n^2)),
  knot_y      = factor(rep(idx[, 2], each = n^2)),
  basis_fun   = factor(rep(1:(df^2), each = n^2)),
  x           = covariates[, 1],
  y           = covariates[, 2]
)

dat %>%
  ggplot(aes(x = x, y = spline, color = basis_fun, group = basis_fun)) +
  facet_wrap(~ y, ncol = 5) +
  theme(legend.position = "none") +
  geom_line()

## simulate some data

# The coefficients
Q       <- spline_interactions_penalty(df, phi = 0.99)
beta <- drop(rmvn(1, rep(0, nrow(Q)), solve(Q)))
dat_beta <- data.frame(beta = beta, knot_x = idx[, 1], knot_y = idx[, 2])
dat_beta %>%
  ggplot(aes(x = knot_x, y = knot_y, fill = beta)) +
  scale_fill_viridis_c() +
  geom_raster()

# The mean surface
sigma <- 0.25
mu <- drop(X_bs %*% beta)
y <- mu + rnorm(n^2, 0, sigma)
dat <- data.frame(
  mu = mu,
  y  = y,
  x1 = covariates[, 1],
  x2 = covariates[, 2]
)
dat %>%
  ggplot(aes(x = x1, y = x2, fill = mu)) +
  scale_fill_viridis_c() +
  geom_raster()

dat %>%
  ggplot(aes(x = x1, y = x2, fill = y)) +
  scale_fill_viridis_c() +
  geom_raster()

