library(tidyverse)
library(patchwork)
library(here)
library(rstan)
library(splines)
library(mvnfast)
library(bayesplot)
library(mvnfast)
if (!require(devtools)) {
  install.packages("devtools")
}
if (!require(treeRingSplines)) {
  devtools::install_github("jtipton25/treeRingSplines")
}

options(mc.cores = parallel::detectCores())

# setup directories

if (!dir.exists(here::here("images")))
  dir.create(here::here("images"))
if (!dir.exists(here::here("images", "spline-interaction")))
  dir.create(here::here("images", "spline-interaction"))


# B-spline model with interactions ---------------------------------------------

## simulate some example data with non-linear predictors using the tree ring model

set.seed(99)
n_tree_per_plot <- 10
n_plot <- 20
n_tree <- n_tree_per_plot * n_plot
## one linear predictor
K <- 1
## three non-linear predictors - age, climate normal, yearly variation
q <- 3
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
X <- rep(rnorm(n_plot), each = n_per_tree * n_tree_per_plot)
beta <- rnorm(1)

# construct the spline basis
df <- 5

## a n by q matrix of coefficients for the basis expansion
Z <- cbind(
  # tree-specific predictor over age
  rep(seq(0, 1, length.out = n_per_tree), times = n_tree),
  # plot-level predictor (like climate normal)
  rep(rnorm(n_plot), each = n_per_tree * n_tree_per_plot),
  # plot-level predictor (like climate annual variation)
  rep(rnorm(n_plot * n_per_tree), each = n_tree_per_plot)
)

## build the univariate splines from a matrix of covariates
X_bs  <- build_spline(Z, df)

sigma_beta_bs <- c(0.5, 0.01, 0.1)

## figure out how to do this more programmatically (loop/function)
beta_bs <- matrix(0, 3, df)
for (j in 1:3) {
  beta_bs[j, 1] <- rnorm(1, 0, sigma_beta_bs[j])
}

## Autoregressive prior on beta's to induce smoothing
for (i in 2:df) {
  for (j in 1:3) {
    beta_bs[j, i] <- rnorm(1, beta_bs[, i-1], sigma_beta_bs[j])
  }
}

## plot the basis functions
data.frame(
  x           = Z[, 1],
  y           = c(X_bs[1, , ]),
  observation = rep(1:nrow(X_bs[1, , ]), times = ncol(X_bs[1, , ])),
  basis_fun   = rep(1:ncol(X_bs[1, , ]), each = nrow(X_bs[1, , ]))
) %>%
  ggplot(aes(x = x, y = y, group = basis_fun, color = basis_fun)) +
  geom_line()

# plot the regression coefficients as a function of the knot
data.frame(
  beta = c(beta_bs),
  par = factor(rep(1:3, times = df)),
  knot = rep(1:df, each = 3)
) %>%
  ggplot(aes(x = knot, y = beta, group = par, color = par)) +
  geom_line() +
  ggtitle("Simulated betas")


# plot the marginal simulated spline effects
data.frame(
  effect = c(build_spline_effects(X_bs, beta_bs)),
  par        = factor(rep(1:3, each = n)),
  X          = c(Z)
) %>%
  ggplot(aes(x = X, y = effect, group = par)) +
  geom_line() +
  facet_wrap(~ par)

## Build the spline interactions
X_bs_int  <- build_spline_interactions(Z, df)
## spline interaction smoothness penalty (2D)
## this is set at phi = 0.99 for simulation, for fitting in stan, we set phi = 1
Q <- make_Q(df, phi = 0.99, prec_model = "CAR")


beta_bs_int <- matrix(0, nrow(X_bs_int), df^2)
sigma_beta_bs_int <- c(0.05, 0.01, 0.025)
for (j in 1:nrow(X_bs_int)) {
  beta_bs_int[j, ] <- drop(rmvn(1, rep(0, df^2), sigma_beta_bs_int[j] * solve(Q)))
}




sigma <- 0.1
## the drop() makes y into a vector
y <- drop(beta0 + as.matrix(X) %*% beta + rowSums(build_spline_effects(X_bs, beta_bs)) + rowSums(build_spline_effects(X_bs_int, beta_bs_int)) + rnorm(n, 0, sigma))



dat <- data.frame(
  y = y,
  X = X,
  X1 = Z[, 1],
  X2 = Z[, 2],
  X3 = Z[, 3],
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
  ggplot(aes(x = X, y = y, color = plot)) +
  geom_point(alpha = 0.25) +
  ggtitle("tree increment as a function of linear variable") +
  scale_colour_viridis_d(end = 0.8) +
  theme(legend.position = "none") +
  xlab("linear variable")

p3 <- dat %>%
  ggplot(aes(x = X2, y = y, color = plot)) +
  geom_point(alpha = 0.25) +
  ggtitle("tree increment as a function of climate normal") +
  scale_colour_viridis_d(end = 0.8) +
  theme(legend.position = "none") +
  xlab("climate normal")


p4 <- dat %>%
  ggplot(aes(x = X3, y = y, color = plot)) +
  geom_point(alpha = 0.125) +
  ggtitle("tree increment as a function of annual climate") +
  theme(legend.position = "none") +
  scale_colour_viridis_d(end = 0.8) +
  xlab("annual climate")

(p1 + p2) / (p3 + p4)

## Plot the marginal and interaction effects

n_grid <- 40
Z_grid <- expand.grid(
  seq(min(Z[, 1]), max(Z[, 1]), length.out = n_grid),
  seq(min(Z[, 2]), max(Z[, 2]), length.out = n_grid),
  seq(min(Z[, 3]), max(Z[, 3]), length.out = n_grid)
)
Z_bs_grid <- build_spline(as.matrix(Z_grid), df)
effects_grid <- build_spline_effects(Z_bs_grid, beta_bs)

Z_bs_int_grid <- build_spline_interactions(as.matrix(Z_grid), df)
effects_int_grid <- build_spline_effects(Z_bs_int_grid, beta_bs_int)

dat_marginal <- data.frame(
  x1 = Z_grid[, 1],
  x2 = Z_grid[, 2],
  x3 = Z_grid[, 3],
  effects1 = effects_grid[, 1],
  effects2 = effects_grid[, 2],
  effects3 = effects_grid[, 3]
  )

p1 <- dat_marginal %>%
  ggplot(aes(x = x1, y = effects1)) +
  geom_line() +
  ggtitle("marginal effect")
p2 <- dat_marginal %>%
  ggplot(aes(x = x2, y = effects2)) +
  geom_line() +
  ggtitle("marginal effect")
p3 <- dat_marginal %>%
  ggplot(aes(x = x3, y = effects3)) +
  geom_line() +
  ggtitle("marginal effect")

dat_int <- data.frame(
  x1 = Z_grid[, 1],
  x2 = Z_grid[, 2],
  x3 = Z_grid[, 3],
  effects1 = effects_int_grid[, 1],
  effects2 = effects_int_grid[, 2],
  effects3 = effects_int_grid[, 3]
)

p4 <- dat_int %>%
  ggplot(aes(x = x1, y = x2, fill = effects1)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("interaction")
p5 <- dat_int %>%
  ggplot(aes(x = x1, y = x3, fill = effects2)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("interaction")
p6 <- dat_int %>%
  ggplot(aes(x = x2, y = x3, fill = effects3)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("interaction")

dat_full <- data.frame(
  x1 = Z_grid[, 1],
  x2 = Z_grid[, 2],
  x3 = Z_grid[, 3],
  effects1 = effects_grid[, 1] + effects_grid[, 2] + effects_int_grid[, 1],
  effects2 = effects_grid[, 1] + effects_grid[, 3] + effects_int_grid[, 2],
  effects3 = effects_grid[, 2] + effects_grid[, 3] + effects_int_grid[, 3]
)

p7 <- dat_full %>%
  ggplot(aes(x = x1, y = x2, fill = effects1)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("full effect")
p8 <- dat_full %>%
  ggplot(aes(x = x1, y = x3, fill = effects2)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("full effect")
p9 <- dat_full %>%
  ggplot(aes(x = x2, y = x3, fill = effects3)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("full effect")

(p1 + p2 + p3) / (p4 + p5 + p6) / (p7 + p8 + p9)

##----------------------------------------------------------------------------##
## Fit the b-spline model with interactions in stan
##----------------------------------------------------------------------------##

## for now, predict at the observed sites only
X_pred <- X
X_bs_pred <- X_bs
X_bs_int_pred <- X_bs_int
tree_idx_pred <- tree_idx

## spline interaction smoothness penalty (2D)
## this is set at phi = 0.99 for simulation, for fitting in stan, we set phi = 1
Q <- make_Q(df, phi = 0.99, prec_model = "CAR")

## Needs to fit with more samples
if (file.exists(here("results", "spline-interaction-example.RDS"))) {
  fit_grow <- readRDS(here("results", "spline-interaction-example.RDS"))
} else {
  fit_grow <- lm_splines_interaction(
    y                = y,
    X                = as.matrix(X),
    X_bs             = X_bs,
    X_bs_int         = X_bs_int,
    Q                = Q,
    n_plot           = n_plot,
    n_tree           = n_tree,
    plot_by_tree_idx = plot_by_tree_idx,
    tree_idx         = tree_idx,
    X_pred           = as.matrix(X_pred),
    X_bs_pred        = X_bs_pred,
    X_bs_int_pred    = X_bs_int_pred,
    tree_idx_pred    = tree_idx_pred,
    iter = 2000,
    warmup = 1000,
    chains = 4,
    control = list(max_treedepth = 15, adapt_delta = 0.99)
  )

  saveRDS(fit_grow, file = here("results" ,"spline-interaction-example.RDS"))
}


pars <- rstan::extract(fit_grow)


# check sampling diagnostics
check_hmc_diagnostics(fit_grow)

# check trace plots
p1 <- mcmc_trace(fit_grow, pars = "lp__")
p2 <- mcmc_trace(fit_grow, pars = "mu_beta0")
p3 <- mcmc_trace(fit_grow, pars = vars(param_range("beta", 1:3)))
p4 <- mcmc_trace(fit_grow, regex_pars = "s_|sigma")

## examine MCMC algorithm
mcmc_pairs(
  fit_grow,
  np = nuts_params(fit_grow),
  pars = c("lp__", "sigma_beta[1]", "sigma_beta[2]", "sigma_beta[3]")
)

if (!file.exists(here::here("images", "spline", "spline-trace-others.png"))) {
  ggsave(
    file = here::here("images", "spline", "spline-trace-others.png"),
    width = 16,
    height = 9,
    (p1 + p2 + p3) / p4
  )
}


# check trace plots for tree level intercept
for (j in 1:5) {
  if (!file.exists(here::here("images", "spline", paste0("spline-trace-beta0_t-", j, ".png")))) {
    ggsave(
      file = here::here("images", "spline", paste0("spline-trace-beta0_t-", j, ".png")),
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
if (!file.exists(here::here("images", "spline", "spline-trace-beta0_p.png"))) {
  ggsave(
    file = here::here("images", "spline", "spline-trace-beta0_p.png"),
    width = 16,
    height = 9,
    mcmc_trace(fit_grow,
               pars = vars(param_range("beta0_p", 1:20)),
               facet_args = list(ncol = 2, nrow = 10)) +
      theme_bw(base_size = 14)
  )
}


# check trace plots plot level intercepts
if (!file.exists(here::here("images", "spline", "spline-trace-betas.png"))) {
  p1 <- mcmc_trace(fit_grow, pars = "beta[1]") + #vars(param_range("beta", 1))) +
    theme_bw(base_size = 14)
  p2 <- mcmc_trace(fit_grow, pars = vars(param_glue("beta_bs[{level1},{level2}]",
                                                    level1 = 1:3, level2 = 1:5)),
                   facet_args = list(nrow = 5, ncol = 3)) +
    theme_bw(base_size = 14)

  ggsave(
    file = here::here("images", "spline", "spline-trace-betas.png"),
    width = 16,
    height = 9,
    p1 / p2 + plot_layout(heights = c(1, 4))
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
data.frame(
  beta      = pars$beta
) %>%
  ggplot(aes(x = "", y = beta)) +
  geom_boxplot() +
  geom_point(position = "jitter", alpha = 0.1) +
  geom_point(data = data.frame(beta = beta), color = "red") +
  ggtitle("slope estimate")

# b-spline parameters
beta_bs_post <- pars$beta_bs
dimnames(beta_bs_post) <- list(
  iteration = 1:dim(beta_bs_post)[1],
  parameter = 1:dim(beta_bs_post)[2],
  knot      = 1:dim(beta_bs_post)[3]
)
as.data.frame.table(beta_bs_post, responseName = "beta") %>%
  mutate(
    iteration = factor(iteration),
    parameter = factor(parameter),
    knot      = factor(knot)
  ) %>%
  group_by(parameter, knot) %>%
  summarize(
    estimate   = mean(beta),
    lower_beta = quantile(beta, prob = 0.025),
    upper_beta = quantile(beta, prob = 0.975)
  ) %>%
  left_join(
    data.frame(
      truth     = c(beta_bs),
      parameter = factor(rep(1:3, times = df)),
      knot      = factor(rep(1:df, each = 3))
    )
  ) %>%
  ggplot(aes(x = truth, y = estimate, color = knot)) +
  geom_point(alpha = 0.75) +
  geom_errorbar(aes(ymin = lower_beta, ymax = upper_beta), alpha = 0.5, width = 0.0) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_wrap(~ parameter) +
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



data.frame(
  sigma_beta = c(pars$sigma_beta),
  sigma_true = rep(sigma_beta_bs, each = nrow(pars$sigma_beta)),
  par        = factor(rep(1:length(sigma_beta_bs), each = nrow(pars$sigma_beta)))
) %>%
  ggplot(aes(x = par, y = sigma_beta)) +
  geom_violin() +
  geom_point(alpha = 0.1, position = "jitter") +
  geom_point(aes(x = par, y = sigma_true), color = "red")  +
  ggtitle("spline penalty variance estimate")



## Fitted vs. estimated functional responses

dat_X <- data.frame(
  X           = c(Z[, 1], Z[, 2], Z[, 3]),
  observation = 1:length(dat$y),
  parameter   = factor(rep(paste("var", 1:3, sep="-"), each = length(dat$y)))
)

# dat_data <- data.frame(
#   y = dat$y,
#
# )

effects <- array(
  0,
  dim      = c(nrow(pars$beta), length(y), 3),
  dimnames = list(
    iteration   = 1:nrow(pars$beta),
    observation = 1:length(y),
    parameter   = paste("var", 1:3, sep="-")
  )
)

for (k in 1:nrow(pars$beta)) {
  effects[k, , 1] <- X_bs[1, , ] %*% pars$beta_bs[k, 1, ]
  effects[k, , 2] <- X_bs[2, , ] %*% pars$beta_bs[k, 2, ]
  effects[k, , 3] <- X_bs[3, , ] %*% pars$beta_bs[k, 3, ]
}

dat_effects <- as.data.frame.table(effects, responseName = "effect")
dat_effects$iteration <- as.numeric(dat_effects$iteration)
dat_effects$observation <- as.numeric(dat_effects$observation)
dat_effects <- dat_effects %>%
  left_join(dat_X)

dat_truth <- data.frame(
  effect      = c(
    X_bs[1, , ] %*% beta_bs[1, ],
    X_bs[2, , ] %*% beta_bs[2, ],
    X_bs[3, , ] %*% beta_bs[3, ]
  ),
  X           = c(Z[, 1], Z[, 2], Z[, 3]),
  observation = 1:length(dat$y),
  parameter   = factor(rep(paste("var", 1:3, sep="-"), each = length(dat$y)))

)


if (!file.exists(here::here("images", "spline", "spline-effects.png"))) {
  ggsave(
    file = here::here("images", "spline", "spline-effects.png"),
    width = 16,
    height = 9,
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
  )
}

# B-spline model with spline interactions --------------------------------------









# ## make sure the igraph and Matrix packages are installed
#
# make_Q <- function(n, phi) {
#   # if (n <= 2)
#   #     stop("n must be an integer larger or equal to 2")
#   if (phi < -1 || phi > 1)
#     stop("phi must be between -1 and 1")
#   W <- igraph::as_adjacency_matrix(
#     igraph::make_lattice(
#       length = n,
#       dim    = 1
#     ),
#     sparse = TRUE
#   )
#   D <- Matrix::Diagonal(x = Matrix::colSums(W))
#   Q <- D - phi * W
#   return(Q)
# }
#
# Q_beta <- sapply(1:length(df), function(i) as.matrix(make_Q(df[i], 0.9)), simplify = FALSE)
