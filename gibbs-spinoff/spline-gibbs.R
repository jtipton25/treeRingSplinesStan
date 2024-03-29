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
  devtools::install_github("jtipton25/treeRingSplineStan")
}

options(mc.cores = parallel::detectCores())

# setup directories

if (!dir.exists(here::here("images")))
  dir.create(here::here("images"))
if (!dir.exists(here::here("images", "spline-gibbs")))
  dir.create(here::here("images", "spline-gibbs"))


# B-spline model --------------------------------------------------------------

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
mu_beta0   <- rnorm(1, 0, 0.025)
beta0_tree <- c(0, rnorm(n_tree - 1, 0, 0.1))
beta0_plot <- c(0, rnorm(n_plot - 1, 0, 0.05))
beta0 <- mu_beta0 + beta0_tree[tree_idx] + beta0_plot[plot_idx]

# design matrix representation
beta_intercept <- c(mu_beta0, beta0_plot[-1], beta0_tree[-1])
X_intercept <- model.matrix(~ factor(plot_idx) + factor(tree_idx))
plot(X_intercept %*% beta_intercept, beta0)
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
Z_plot <- rnorm(n_plot * n_per_tree)
site_age_idx <- c(
  rep(1:100, 10), rep(101:200, 10), rep(201:300, 10),
  rep(301:400, 10), rep(401:500, 10), rep(501:600, 10),
  rep(601:700, 10), rep(701:800, 10), rep(801:900, 10),
  rep(901:1000, 10), rep(1001:1100, 10), rep(1101:1200, 10),
  rep(1201:1300, 10), rep(1301:1400, 10), rep(1401:1500, 10),
  rep(1501:1600, 10), rep(1601:1700, 10), rep(1701:1800, 10),
  rep(1801:1900, 10), rep(1901:2000, 10)
)

Z <- cbind(
  # tree-specific predictor over age
  rep(seq(0, 1, length.out = n_per_tree), times = n_tree),
  # plot-level predictor (like climate normal)
  rep(rnorm(n_plot), each = n_per_tree * n_tree_per_plot),
  # plot-level predictor (like climate annual variation)
  Z_plot[site_age_idx]
)

## build the univariate splines from a matrix of covariates
X_bs  <- build_spline(Z, df)
Xbs   <- cbind(X_bs[1, , ], X_bs[2, , ], X_bs[3, , ])

sigma_beta_bs <- c(0.5, 0.01, 0.1)

## figure out how to do this more programmatically (loop/function)
beta_bs <- matrix(0, 3, df)
beta_bs[, 1] <- rnorm(3, 0, sigma_beta_bs)

## Autoregressive prior on beta's to induce smoothing
for (i in 2:df) {
  beta_bs[, i] <- rnorm(3, beta_bs[, i-1], sigma_beta_bs)
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


# plot the simulated spline effects
data.frame(
  effect = c(build_spline_effects(X_bs, beta_bs)),
  par        = factor(rep(1:3, each = n)),
  X          = c(Z)
) %>%
  ggplot(aes(x = X, y = effect, group = par)) +
  geom_line() +
  facet_wrap(~ par)



sigma <- 0.1
## the drop() makes y into a vector
y <- drop(beta0 + as.matrix(X) %*% beta + rowSums(build_spline_effects(X_bs, beta_bs)) + rnorm(n, 0, sigma))



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

##----------------------------------------------------------------------------##
## Fit the b-spline model in stan
##----------------------------------------------------------------------------##

## for now, predict at the observed sites only
X_pred <- X
X_bs_pred <- X_bs
tree_idx_pred <- tree_idx

params <- list(
  n_mcmc    = 5000,
  n_adapt   = 5000,
  n_message = 50,
  n_thin    = 10
)


## Needs to fit with more samples
if (file.exists(here("results", "spline-example-gibbs.RDS"))) {
  fit_grow <- readRDS(here("results", "spline-example-gibbs.RDS"))
} else {
  fit_grow <- mcmc_splines(
    y        = y,
    X        = as.matrix(X),
    Xbs      = Xbs,
    plot_idx = plot_idx,
    tree_idx = tree_idx,
    params   = params,
    verbose  = FALSE)

  saveRDS(fit_grow, file = here("results" ,"spline-example-gibbs.RDS"))
}


# check trace plots of log probability
data.frame(
  lp__        = c(fit_grow$lp__),
  observation = rep(1:ncol(fit_grow$lp__), each = nrow(fit_grow$lp__)),
  iteration   = rep(1:nrow(fit_grow$lp__), time = ncol(fit_grow$lp__))
) %>%
  sample_n_of(10, observation) %>%
  ggplot(aes(x = iteration, y = lp__, group = observation, color = observation)) +
  geom_line(alpha = 0.25) +
  theme(legend.position = "none")

# trace plots of intercepts
data.frame(
  beta_int    = c(fit_grow$beta_intercept),
  parameter   = rep(1:ncol(fit_grow$beta_intercept), each = nrow(fit_grow$beta_intercept)),
  iteration   = rep(1:nrow(fit_grow$beta_intercept), time = ncol(fit_grow$beta_intercept))
) %>%
  sample_n_of(10, parameter) %>%
  ggplot(aes(x = iteration, y = beta_int, group = parameter, color = parameter)) +
  geom_line(alpha = 0.25) +
  theme(legend.position = "none")

# trace plots of slopes
data.frame(
  beta        = c(fit_grow$beta),
  parameter   = rep(1:ncol(fit_grow$beta), each = nrow(fit_grow$beta)),
  iteration   = rep(1:nrow(fit_grow$beta), time = ncol(fit_grow$beta))
) %>%
  ggplot(aes(x = iteration, y = beta, group = parameter, color = parameter)) +
  geom_line(alpha = 0.25) +
  theme(legend.position = "none")

# trace plots of spline parameters
data.frame(
  beta        = c(fit_grow$beta_bs),
  parameter   = factor(rep(1:ncol(Z), each = nrow(fit_grow$beta) * df)),
  knot        = factor(rep(1:df, each = nrow(fit_grow$beta))),
  iteration   = rep(1:nrow(fit_grow$beta), time = ncol(fit_grow$beta))
) %>%
  ggplot(aes(x = iteration, y = beta, group = knot, color = knot)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ parameter) +
  theme(legend.position = "none")

# plot of intercept parameters
data.frame(
  intercept_fit = apply(fit_grow$beta_intercept, 2, mean),
  intercept_truth = beta_intercept,
  parameter = c("overall", rep("plot", n_plot - 1), rep("tree", n_tree - 1))
) %>%
  ggplot(aes(x = intercept_fit, y = intercept_truth, color = parameter)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red")

layout(matrix(1:4, 2, 2))
hist(sqrt(fit_grow$sigma2))
abline(v = sigma, col = "red", lwd = 3)
hist(fit_grow$s2_beta_p)
abline(v = 0.05, col = "red", lwd = 3)
hist(fit_grow$s2_beta_t)
abline(v = 0.01, col = "red", lwd = 3)
plot(apply(fit_grow$beta_bs, 2, mean), c(t(beta_bs)))
abline(0, 1, col = "red")

# pseudo-replication
data.frame(
  climate_normal = X,
  plot           = plot_idx
) %>%
  ggplot(aes(x = plot_idx, y = climate_normal)) +
  geom_point(position = "jitter")



## Joint random effects and climate normal effects are identifiable
##  suggest dropping plot level random effects if you want to interpret
## the climate normals!
data.frame(
  estimated_intercept = c(apply(fit_grow$beta_intercept, 2, mean)[1] * rep(1, n_plot) +
                            c(0, apply(fit_grow$beta_intercept, 2, mean)[2:n_plot])
                            )[plot_idx] +
    X * mean(fit_grow$beta) +
    X_bs[2, , ] %*% apply(fit_grow$beta_bs, 2, mean)[df + 1:df],
  sim_intercept = c(mu_beta0 * rep(1, n_plot) + beta0_plot)[plot_idx] +
    X * beta + X_bs[2, , ] %*% beta_bs[2, ]) %>%
  ggplot(aes(x = estimated_intercept, y = sim_intercept)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red")

# co-linearity between plot-level intercepts and climate normals
# shows up again
# plot(
#   apply(fit_grow$beta_intercept, 2, mean)[1:n_plot][plot_idx] + X_bs[2, , ] %*% apply(fit_grow$beta_bs, 2, mean)[df + 1:df],
#   beta0_plot[plot_idx] + X_bs[2, , ] %*% beta_bs[2, ])
# abline(0, 1, col = "red")






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

if (!file.exists(here::here("images", "spline-gibbs", "spline-trace-others.png"))) {
  ggsave(
    file = here::here("images", "spline-gibbs", "spline-trace-others.png"),
    width = 16,
    height = 9,
    (p1 + p2 + p3) / p4
  )
}


# check trace plots for tree level intercept
for (j in 1:5) {
  if (!file.exists(here::here("images", "spline-gibbs", paste0("spline-trace-beta0_t-", j, ".png")))) {
    ggsave(
      file = here::here("images", "spline-gibbs", paste0("spline-trace-beta0_t-", j, ".png")),
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
if (!file.exists(here::here("images", "spline-gibbs", "spline-trace-beta0_p.png"))) {
  ggsave(
    file = here::here("images", "spline-gibbs", "spline-trace-beta0_p.png"),
    width = 16,
    height = 9,
    mcmc_trace(fit_grow,
               pars = vars(param_range("beta0_p", 1:20)),
               facet_args = list(ncol = 2, nrow = 10)) +
      theme_bw(base_size = 14)
  )
}


# check trace plots plot level intercepts
if (!file.exists(here::here("images", "spline-gibbs", "spline-trace-betas.png"))) {
  p1 <- mcmc_trace(fit_grow, pars = "beta[1]") + #vars(param_range("beta", 1))) +
    theme_bw(base_size = 14)
  p2 <- mcmc_trace(fit_grow, pars = vars(param_glue("beta_bs[{level1},{level2}]",
                                                    level1 = 1:3, level2 = 1:5)),
                   facet_args = list(nrow = 5, ncol = 3)) +
    theme_bw(base_size = 14)

  ggsave(
    file = here::here("images", "spline-gibbs", "spline-trace-betas.png"),
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


if (!file.exists(here::here("images", "spline-gibbs", "spline-effects.png"))) {
  ggsave(
    file = here::here("images", "spline-gibbs", "spline-effects.png"),
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
