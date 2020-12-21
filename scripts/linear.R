library(tidyverse)
library(patchwork)
library(here)
library(rstan)
library(splines)
library(mvnfast)
if (!require(devtools)) {
  install.packages("devtools")
}
if (!require(treeRingSplines)) {
  devtools::install_github("jtipton25/treeRingSplines")
}

options(mc.cores = parallel::detectCores())

# B-spline model --------------------------------------------------------------

## simulate some example data with non-linear predictors using the tree ring model

set.seed(99)
n_tree_per_plot <- 10
n_plot <- 20
n_tree <- n_tree_per_plot * n_plot
## three linear predictors - age, climate normal, yearly variation
K <- 3
# observations per tree
n_per_tree <- 100
## number of observations
n  <- n_tree * n_per_tree
n_pred <- n_tree * n_per_tree

# plot index
plot_idx <- rep(1:n_plot, each = n_tree_per_plot * n_per_tree)
tree_idx <- rep(1:n_tree, each = n_per_tree)
plot_by_tree_idx <- rep(1:n_plot, each = n_tree_per_plot)

# intercept terms
beta0_tree <- rnorm(n_tree, 0, 0.1)
beta0_plot <- rnorm(n_plot, 0, 0.05)
beta0 <- beta0_tree[tree_idx] + beta0_plot[plot_idx]

# plot the intercepts
data.frame(
  beta0      = beta0,
  beta0_plot = beta0_plot[plot_idx],
  tree        = factor(tree_idx),
  plot        = factor(plot_idx)
) %>%
  ggplot(aes(x = plot, y = beta0)) +
  geom_boxplot() +
  geom_point(alpha = 0.005) +   ## low number due to replication of points
  geom_point(aes(x = plot, y = beta0_plot), color = "red") +
  ggtitle("Tree and plot level intercepts")

# construct an example linear term (one per plot)
X <- cbind(
  # tree-specific predictor over age
  rep(seq(0, 1, length.out = n_per_tree), times = n_tree),
  # plot-level predictor (like climate normal)
  rep(rnorm(n_plot), each = n_per_tree * n_tree_per_plot),
  # plot-level predictor (like climate annual variation)
  rep(rnorm(n_plot * n_per_tree), each = n_tree_per_plot)
)

beta <- rnorm(ncol(X))


# plot the simulated spline effects
data.frame(
  effect = c(X %*% beta),
  par        = factor(rep(1:3, each = n)),
  X          = c(X)
) %>%
  ggplot(aes(x = X, y = effect, group = par)) +
  geom_line() +
  facet_wrap(~ par)

## measurement error
sigma <- 0.1
## the drop() makes y into a vector
y <- drop(beta0 + X %*% beta + rnorm(n, 0, sigma))


## plot the simulated data
dat <- data.frame(
  y = y,
  X1 = X[, 1],
  X2 = X[, 2],
  X3 = X[, 3],
  tree = factor(tree_idx),
  plot = factor(plot_idx),
  time = 1:100
)

p1 <- dat %>%
  ggplot(aes(x = time, y = y, group = tree)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ plot) +
  ggtitle("example simulated data")

p2 <- dat %>%
  ggplot(aes(x = X2, y = y, color = plot)) +
  geom_point(alpha = 0.25) +
  ggtitle("tree increment as a function of climate normal") +
  scale_colour_viridis_d(end = 0.8) +
  theme(legend.position = "none") +
  xlab("climate normal")

p3 <- dat %>%
  ggplot(aes(x = X3, y = y, color = plot)) +
  geom_point(alpha = 0.125) +
  ggtitle("tree increment as a function of annual climate") +
  theme(legend.position = "none") +
  scale_colour_viridis_d(end = 0.8) +
  xlab("annual climate")

p1 / (p2 + p3)

##----------------------------------------------------------------------------##
## Fit the linear model in stan
##----------------------------------------------------------------------------##

## for now, predict at the observed sites only
X_pred <- X
tree_idx_pred <- tree_idx

## Needs to fit with more samples
if (file.exists(here("results", "linear-example.RDS"))) {
  fit_grow <- readRDS(here("results", "linear-example.RDS"))
} else {
  fit_grow <- lm_linear(
    y                = y,
    X                = X,
    n_plot           = n_plot,
    n_tree           = n_tree,
    plot_by_tree_idx = plot_by_tree_idx,
    tree_idx         = tree_idx,
    X_pred           = X_pred,
    tree_idx_pred    = tree_idx_pred,
    iter = 1000,
    warmup = 500,
    chains = 4,
    control = list(max_treedepth = 15)
  )

  saveRDS(fit_grow, file = here("results" ,"linear-example.RDS"))
}


pars <- rstan::extract(fit_grow)

## plot the estimated vs. fitted intercepts
data.frame(
  beta0     = c(pars$beta0_t),
  iteration = rep(1:nrow(pars$beta0_t), times = n_tree),
  tree      = factor(rep(1:n_tree, each = nrow(pars$beta0_t))),
  plot      = factor(rep(plot_by_tree_idx, each = nrow(pars$beta0_t)))
) %>%
  group_by(tree, plot) %>%
  summarize(
    estimate    = mean(beta0),
    lower_beta0 = quantile(beta0, prob = 0.025),
    upper_beta0 = quantile(beta0, prob = 0.975)
  ) %>%
  left_join(
    data.frame(
      truth = beta0[seq(1, length(beta0), by = n_per_tree)],
      tree  = factor(1:n_tree),
      plot  = factor(plot_by_tree_idx)
    )
  ) %>%
  ggplot(aes(x = truth, y = estimate)) +
  geom_point(alpha = 0.75) +
  geom_errorbar(aes(ymin = lower_beta0, ymax = upper_beta0), alpha = 0.5, width = 0) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("tree level intercept estimates\n(including plot level effects)")

## plot the estimated vs. fitted regression coefficients
beta_post <- pars$beta
dimnames(beta_post) <- list(
  iteration = 1:dim(beta_post)[1],
  parameter = 1:dim(beta_post)[2]
)
as.data.frame.table(beta_post, responseName = "beta") %>%
  mutate(
    iteration = factor(iteration),
    parameter = factor(parameter)
  ) %>%
  group_by(parameter) %>%
  summarize(
    estimate   = mean(beta),
    lower_beta = quantile(beta, prob = 0.025),
    upper_beta = quantile(beta, prob = 0.975)
  ) %>%
  left_join(
    data.frame(
      truth     = c(beta),
      parameter = factor(1:3)
    )
  ) %>%
  ggplot(aes(x = truth, y = estimate)) +
  geom_point(alpha = 0.75) +
  geom_errorbar(aes(ymin = lower_beta, ymax = upper_beta), alpha = 0.5, width = 0.0) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("slope estimates")


data.frame(
  sigma = pars$sigma_y,
  sigma_true = sigma
) %>%
  ggplot(aes(x = "", y = sigma)) +
  geom_violin() +
  geom_point(alpha = 0.1, position = "jitter") +
  geom_point(aes(x = "", y = sigma_true), color = "red")  +
  ggtitle("observation variance estimate")




## Fitted vs. estimated functional responses

dat_X <- data.frame(
  X           = c(X),
  observation = factor(1:length(dat$y)),
  parameter   = factor(rep(paste("var", 1:3, sep="-"), each = length(dat$y)))
)
effects <- array(0, dim = c(ncol(X), nrow(X), nrow(pars$beta)),
                 dimnames = list(
                   parameter   = paste("var", 1:ncol(X), sep = '-'),
                   observation = 1:nrow(X),
                   iteration   = 1:nrow(pars$beta)
                 )
)
effects[1, , ] <- X[, 1, drop = FALSE] %*% t(pars$beta[, 1])
effects[2, , ] <- X[, 2, drop = FALSE] %*% t(pars$beta[, 2])
effects[3, , ] <- X[, 3, drop = FALSE] %*% t(pars$beta[, 3])

dat_effects <- as.data.frame.table(effects, responseName = "effect") %>%
  mutate(
    iteration = factor(iteration),
    parameter = factor(parameter)
  ) %>%
  left_join(dat_X)

dat_truth <- data.frame(
  effect      = c(X[, 1] * beta[1], X[, 2] * beta[2], X[, 3] * beta[3]),
  X           = c(X[, 1], X[, 2], X[, 3]),
  observation = 1:length(dat$y),
  parameter   = factor(rep(paste("var", 1:3, sep="-"), each = n))
)

dat_effects %>%
  group_by(observation, parameter, X) %>%
  summarize(
    effect_mean = mean(effect),
    effect_lower = quantile(effect, prob = 0.025),
    effect_upper = quantile(effect, prob = 0.975)
  ) %>%
  ggplot(aes(x = X, y = effect_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = effect_lower, ymax = effect_upper), fill = "grey", alpha = 0.5) +
  geom_line(data = dat_truth, aes(x = X, y = effect), color = "red") +
  facet_wrap(~ parameter, ncol = 1, scales = "free") +
  ggtitle("estimate in grey, simulated trend in red")



