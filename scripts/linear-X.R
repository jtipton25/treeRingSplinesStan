library(tidyverse)
library(patchwork)
library(here)
library(rstan)
library(splines)
library(mvnfast)
library(bayesplot)
if (!require(devtools)) {
  install.packages("devtools")
}
if (!require(treeRingSplinesStan)) {
  devtools::install_github("jtipton25/treeRingSplinesStan")
}

options(mc.cores = parallel::detectCores())

# setup directories

if (!dir.exists(here::here("images")))
  dir.create(here::here("images"))
if (!dir.exists(here::here("images", "linear-no-p")))
  dir.create(here::here("images", "linear-no-p"))


# linear model --------------------------------------------------------------

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
beta0 <- beta0_tree[tree_idx]

# plot the intercepts
data.frame(
  beta0      = beta0,
  tree        = factor(tree_idx),
  plot        = factor(plot_idx)
) %>%
  ggplot(aes(x = plot, y = beta0)) +
  geom_boxplot() +
  geom_point(alpha = 0.005) +   ## low number due to replication of points
  ggtitle("Tree and plot level intercepts")

# plot-level predictor (like climate annual variation)
X_plot <- rnorm(n_plot * n_per_tree)
site_age_idx <- c(
  rep(1:100, 10), rep(101:200, 10), rep(201:300, 10),
  rep(301:400, 10), rep(401:500, 10), rep(501:600, 10),
  rep(601:700, 10), rep(701:800, 10), rep(801:900, 10),
  rep(901:1000, 10), rep(1001:1100, 10), rep(1101:1200, 10),
  rep(1201:1300, 10), rep(1301:1400, 10), rep(1401:1500, 10),
  rep(1501:1600, 10), rep(1601:1700, 10), rep(1701:1800, 10),
  rep(1801:1900, 10), rep(1901:2000, 10)
)

X <- cbind(
  # tree-specific predictor over age
  rep(seq(0, 1, length.out = n_per_tree), times = n_tree),
  # plot-level predictor (like climate normal)
  rep(rnorm(n_plot), each = n_per_tree * n_tree_per_plot),
  X_plot[site_age_idx]
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
if (file.exists(here("results", "linear-no-p-example.RDS"))) {
  fit_grow <- readRDS(here("results", "linear-no-p-example.RDS"))
} else {
  fit_grow <- lm_linear_no_p(
    y                = y,
    X                = X,
    n_tree           = n_tree,
    tree_idx         = tree_idx,
    X_pred           = X_pred,
    tree_idx_pred    = tree_idx_pred,
    iter = 1000,
    warmup = 500,
    chains = 4,
    control = list(max_treedepth = 15)
  )

  saveRDS(fit_grow, file = here("results" ,"linear-no-p-example.RDS"))
}

pars <- rstan::extract(fit_grow)


print(get_elapsed_time(fit_grow))
print(get_elapsed_time(fit_grow2))

# check sampling diagnostics
check_hmc_diagnostics(fit_grow)

# check trace plots
p1 <- mcmc_trace(fit_grow, pars = "lp__")
p2 <- mcmc_trace(fit_grow, pars = "mu_beta0")
p3 <- mcmc_trace(fit_grow, regex_pars = "s_|sigma")
p4 <- mcmc_trace(fit_grow, pars = vars(param_range("beta", 1:3)))
if (!file.exists(here::here("images", "linear-no-p", "linear-trace-others.png"))) {
  ggsave(
    file = here::here("images", "linear-no-p", "linear-trace-others.png"),
    width = 16,
    height = 9,
    (p1 + p2) / p3 / p4
  )
}


# check trace plots for tree level intercept
for (j in 1:5) {
  if (!file.exists(here::here("images", "linear-no-p", paste0("linear-trace-beta0_t-", j, ".png")))) {
    ggsave(
      file = here::here("images", "linear-no-p", paste0("linear-trace-beta0_t-", j, ".png")),
      width = 16,
      height = 9,
      mcmc_trace(fit_grow,
                 pars = vars(param_range("beta0_t", ((j-1) * 40 + 1):(40 + (j-1) * 40))),
                 facet_args = list(ncol = 4, nrow = 10)) +
        theme_bw(base_size = 14)
    )
  }
}

# check trace plots plot level intercepts
if (!file.exists(here::here("images", "linear-no-p", "linear-trace-beta0_p.png"))) {
  ggsave(
    file = here::here("images", "linear-no-p", "linear-trace-beta0_p.png"),
    width = 16,
    height = 9,
    mcmc_trace(fit_grow,
               pars = vars(param_range("beta0_p", 1:20)),
               facet_args = list(ncol = 2, nrow = 10)) +
      theme_bw(base_size = 14)
  )
}

# check trace plots plot level intercepts
if (!file.exists(here::here("images", "linear-no-p", "linear-trace-beta0_p.png"))) {
  ggsave(
    file = here::here("images", "linear-no-p", "linear-trace-beta0_p.png"),
    width = 16,
    height = 9,
    mcmc_trace(fit_grow,
               pars = vars(param_range("beta0_p", 1:20)),
               facet_args = list(ncol = 2, nrow = 10)) +
      theme_bw(base_size = 14)
  )
}


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
  ggplot(aes(x = truth, y = estimate, color = parameter)) +
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


if (!file.exists(here::here("images", "linear-no-p", "linear-effects.png"))) {

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


  p1 <- dat_effects %>%
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

  ggsave(
    file = here::here("images", "linear-no-p", "linear-effects.png"),
    width = 16,
    height = 9,
    p1 + theme_bw(base_size = 14)
  )
}

# Posterior predictive distributions

if (!file.exists(here::here("images", "linear-no-p", "linear-ppc.png"))) {
  p1 <- ppc_dens_overlay(y, pars$y_rep)
  p2 <- ppc_ecdf_overlay(y, pars$y_rep[sample(nrow(pars$y_rep), 250), ])

  ggsave(
    file = here::here("images", "linear-no-p", "linear-ppc.png"),
    width = 16,
    height = 9,
    (p1 + theme_bw(base_size = 14)) / (p2 + theme_bw(base_size = 14))
  )
}


# explore posterior mean of residuals as a function of covariates
if (!file.exists(here::here("images", "linear-no-p", "linear-ppc-by-covariate.png"))) {
  p1 <- ppc_error_scatter_avg_vs_x(y, pars$y_rep, X[, 1], alpha = 0.1) +
    geom_hline(yintercept = 0, color = "red", lty = 2) +
    theme_bw(base_size = 14)
  p2 <- ppc_error_scatter_avg_vs_x(y, pars$y_rep, X[, 2], alpha = 0.1) +
    geom_hline(yintercept = 0, color = "red", lty = 2) +
    theme_bw(base_size = 14)
  p3 <- ppc_error_scatter_avg_vs_x(y, pars$y_rep, X[, 3], alpha = 0.1) +
    geom_hline(yintercept = 0, color = "red", lty = 2) +
    theme_bw(base_size = 14)

  ggsave(
    file = here::here("images", "linear-no-p", "linear-ppc-by-covariate.png"),
    width = 16,
    height = 9,
    p1 / p2 / p3
  )
}


## LOO
loo_linear <- loo(fit_grow, save_psis = TRUE, cores = 2)
psis_linear <- loo_linear$psis_object
lw <- weights(psis_linear)

# marginal predictive check using LOO probability integral transform
ppc_loo_pit_overlay(y, pars$y_rep, lw = lw)
ppc_loo_pit_qq(y, pars$y_rep, lw = lw)

# loo predictive intervals vs observations
keep_obs <- sample(length(y), 100)
ppc_loo_intervals(y, pars$y_rep, psis_object = psis_linear, subset = keep_obs, order = "median")


plot(loo_linear, main = "PSIS diagnostic plot linear")

