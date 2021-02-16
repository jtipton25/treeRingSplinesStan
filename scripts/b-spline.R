# Things to add
# hierarchical pooling of tree-level random effects 

# fixed double (non-identifiable) intercept mu_beta_0

library(tidyverse)
library(patchwork)
library(here)
library(rstan)
library(splines)
# library(igraph)
# library(Matrix)
library(mvnfast)

options(mc.cores = parallel::detectCores())

# linear model example ---------------------------------------------------------

# simulate some example data with linear predictors using the tree ring model
set.seed(999)
n_tree_per_plot <- 10
n_plot <- 20
ntree <- n_tree_per_plot * n_plot
## three predictors - age, climate normal, yearly variation
K <- 3
# observations per tree
n_per_tree <- 100
## number of observations
nG <- ntree * n_per_tree
n_pred <- ntree * n_per_tree

# plot index
plot_idx <- rep(1:n_plot, each = n_tree_per_plot * n_per_tree)
tree_idx <- rep(1:ntree, each = n_per_tree)
plot_for_tree_idx <- rep(1:n_plot, each = n_tree_per_plot)

# intercept terms
beta0_tree <- rnorm(ntree, 0, 0.1)
beta0_plot <- rnorm(n_plot, 0, 0.05)
beta0 <- beta0_tree[tree_idx] + beta0_plot[plot_idx]

dat_intercepts <- data.frame(
  beta0      = beta0,
  beta0_plot = beta0_plot[plot_idx],
  tree        = factor(tree_idx),
  plot        = factor(plot_idx)
)

dat_intercepts %>%
  ggplot(aes(x = plot, y = beta0)) +
  geom_boxplot() +
  geom_point(alpha = 0.005) +   ## low number due to replication of points
  geom_point(aes(x = plot, y = beta0_plot), color = "red") +
  ggtitle("Tree and plot level intercepts")

# add in tree-specific predictor over age
X1 <- rep(seq(0, 1, length.out = n_per_tree), times = ntree)
beta1 <- 0.1
# add in plot-level predictor (like climate normal)
X2 <- rep(rnorm(n_plot), each = n_per_tree * n_tree_per_plot)
beta2 <- -0.2
# add in plot-level predictor (like climate annual variation)
X3 <- rep(rnorm(n_plot * n_per_tree), each = n_tree_per_plot)
beta3 <- 0.15

sigma <- 0.1
y <- beta0 + X1 * beta1 + X2 * beta2 + X3 * beta3 + rnorm(nG, 0, sigma)

dat <- data.frame(
  y = y,
  X1 = X1, 
  X2 = X2,
  X3 = X3,
  tree = factor(tree_idx),
  plot = factor(plot_idx),
  time = 1:100
)

## Plot linear model simulated data --------------------------------------------

p1 <- dat %>%
  ggplot(aes(x = time, y = y, group = tree)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ plot) +
  ggtitle("example simulated data")

p2 <- dat %>%
  ggplot(aes(x = X1, y = y, group = tree)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ plot) +
  ggtitle("tree increment as a function of age") +
  xlab("age")

p3 <- dat %>%
  ggplot(aes(x = X2, y = y, color = plot)) +
  geom_point(alpha = 0.25) +
  ggtitle("tree increment as a function of climate normal") +
  xlab("climate nomral")


p4 <- dat %>%
  ggplot(aes(x = X3, y = y, color = plot)) +
  geom_point(alpha = 0.25) +
  ggtitle("tree increment as a function of annual climate") +
  xlab("annual climate")



(p1 + p2) / (p3 + p4)

## Fit the linear model in stan ------------------------------------------------


# Note: I remove ppc and log-like for loo for simplicity
pied_dat <- list(K = K, nG = nG, yG = y, xG = cbind(X1, X2, X3), plot = plot_idx, 
                 nplot = n_plot, tree = tree_idx, ntree = ntree, 
                 plotfortree = plot_for_tree_idx)


if (file.exists(here("results", "linear-example.RDS"))) {
  fit_grow <- readRDS(here("results", "linear-example.RDS"))
} else {
  fit_grow <- stan(file = here("b-splines", "pied_grow_linear.stan"),
                   data = pied_dat, 
                   iter = 1000, warmup = 500, chains = 3)
  saveRDS(fit_grow, file = here("results" ,"linear-example.RDS"))
  # Bulk ESS warning -- need to increase number of iterations
}


## Plot linear model results ---------------------------------------------------
pars <- rstan::extract(fit_grow)

## plot the estimated vs. fitted intercepts
data.frame(
  beta0     = c(pars$beta0_t), 
  iteration = rep(1:nrow(pars$beta0_t), times = ntree),
  tree      = factor(rep(1:ntree, each = nrow(pars$beta0_t))),
  plot      = factor(rep(plot_for_tree_idx, each = nrow(pars$beta0_t)))
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
      tree  = factor(1:ntree),
      plot  = factor(plot_for_tree_idx)
    )
  ) %>%
  ggplot(aes(x = truth, y = estimate)) +
  geom_point(alpha = 0.75) +
  geom_errorbar(aes(ymin = lower_beta0, ymax = upper_beta0), alpha = 0.5, width = 0) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("tree level intercept estimates\n(including plot level effects)")

## plot the estimated vs. fitted regression coefficients
data.frame(
  beta      = c(pars$u_beta), 
  iteration = rep(1:nrow(pars$u_beta), times = 3),
  parameter = factor(rep(1:3, each = nrow(pars$u_beta)))
) %>%
  group_by(parameter) %>%
  summarize(
    estimate   = mean(beta), 
    lower_beta = quantile(beta, prob = 0.025),
    upper_beta = quantile(beta, prob = 0.975)
  ) %>%
  left_join(
    data.frame(
      truth      = c(beta1, beta2, beta3),
      parameter  = factor(1:3)
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


# B-spline model --------------------------------------------------------------

## simulate some example data with non-linear predictors using the tree ring model

set.seed(99)
n_tree_per_plot <- 10
n_plot <- 20
ntree <- n_tree_per_plot * n_plot
## three predictors - age, climate normal, yearly variation
K <- 3
# observations per tree
n_per_tree <- 100
## number of observations
n  <- ntree * n_per_tree
n_pred <- ntree * n_per_tree

# plot index
plot_idx <- rep(1:n_plot, each = n_tree_per_plot * n_per_tree)
tree_idx <- rep(1:ntree, each = n_per_tree)
plot_for_tree_idx <- rep(1:n_plot, each = n_tree_per_plot)

# intercept terms
beta0_tree <- rnorm(ntree, 0, 0.1)
beta0_plot <- rnorm(n_plot, 0, 0.05)
beta0 <- beta0_tree[tree_idx] + beta0_plot[plot_idx]

dat_intercepts <- data.frame(
  beta0      = beta0,
  beta0_plot = beta0_plot[plot_idx],
  tree        = factor(tree_idx),
  plot        = factor(plot_idx)
)

dat_intercepts %>%
  ggplot(aes(x = plot, y = beta0)) +
  geom_boxplot() +
  geom_point(alpha = 0.005) +   ## low number due to replication of points
  geom_point(aes(x = plot, y = beta0_plot), color = "red") +
  ggtitle("Tree and plot level intercepts")



# construct penalized spline precision matrix for prior ------------------------
## Let's talk more about this, what it does, and how to use it
## K predictors
df <- 5




# add in tree-specific predictor over age
X1 <- rep(seq(0, 1, length.out = n_per_tree), times = ntree)
# add in plot-level predictor (like climate normal)
X2 <- rep(rnorm(n_plot), each = n_per_tree * n_tree_per_plot)
# add in plot-level predictor (like climate annual variation)
X3 <- rep(rnorm(n_plot * n_per_tree), each = n_tree_per_plot)

## set to a small value for computational stability
X1_bs <- splines::bs(X1, df = df, intercept = FALSE)
X2_bs <- splines::bs(X2, df = df, intercept = FALSE)
X3_bs <- splines::bs(X3, df = df, intercept = FALSE)

sigma_beta <- c(0.5, 0.01, 0.1)

## figure out how to do this more programmatically (loop/function)
beta <- matrix(0, 3, df)
beta[, 1] <- rnorm(3, 0, sigma_beta)

## Autoregressive prior on beta's to induce smoothing
for (i in 2:df) {
  beta[, i] <- rnorm(3, beta[, i-1], sigma_beta)
}

## plot the basis functions
dat <- data.frame(
  x           = X1,
  y           = c(X1_bs),
  observation = rep(1:nrow(X1_bs), times = ncol(X1_bs)),
  basis_fun   = rep(1:ncol(X1_bs), each = nrow(X1_bs))
)

dat %>%
  ggplot(aes(x = x, y = y, group = basis_fun, color = basis_fun)) +
  geom_line()


data.frame(
  beta = c(beta),
  par = factor(rep(1:3, each = df)),
  knot = rep(1:df, times = 3)
) %>%
  ggplot(aes(x = knot, y = beta, group = par, color = par)) +
  geom_line() +
  ggtitle("Simulated betas")

data.frame(
  effect = c(X1_bs %*% beta[1, ], X2_bs %*% beta[2, ], X3_bs %*% beta[3, ]),
  par        = factor(rep(1:3, each = nG)),
  X          = c(X1, X2, X3)
) %>%
  ggplot(aes(x = X, y = effect, group = par)) +
  geom_line() +
  facet_wrap(~ par)



sigma <- 0.1
y <- drop(beta0 + X1_bs %*% beta[1, ]+ X2_bs %*% beta[2, ] +  X3_bs %*% beta[3, ] + rnorm(n, 0, sigma))



dat <- data.frame(
  y = y,
  X1 = X1, 
  X2 = X2,
  X3 = X3,
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
  xlab("climate nomral")


p3 <- dat %>%
  ggplot(aes(x = X3, y = y, color = plot)) +
  geom_point(alpha = 0.125) +
  ggtitle("tree increment as a function of annual climate") +
  theme(legend.position = "none") +
  scale_colour_viridis_d(end = 0.8) +
  xlab("annual climate")

p1 / (p2 + p3)

##----------------------------------------------------------------------------##
## Fit the b-spline model in stan
##----------------------------------------------------------------------------##

x_array <- array(dim = c(3, length(y), 5))
x_array[1, , ] <- X1_bs
x_array[2, , ] <- X2_bs
x_array[3, , ] <- X3_bs

# Note: I remove ppc and log-like for loo for simplicity
pied_dat <- list(p = 3, 
                 K = 5,
                 n = length(y),
                 y = y,
                 x = x_array,
                 plot = plot_idx, 
                 nplot = n_plot, 
                 tree = tree_idx,
                 ntree = ntree, 
                 plotfortree = plot_for_tree_idx
)


## Needs to fit with more samples
if (file.exists(here("results", "spline-example.RDS"))) {
  fit_grow <- readRDS(here("results", "spline-example.RDS"))
} else {
  fit_grow <- stan(file = here("b-splines", "pied_grow_splines.stan"), data = pied_dat, 
                   iter = 1000, warmup = 500, chains = 3, 
                   control = list(max_treedepth = 15)
                   )
  saveRDS(fit_grow, file = here("results" ,"spline-example.RDS"))
}


pars <- rstan::extract(fit_grow)

## plot the estimated vs. fitted intercepts
data.frame(
  beta0     = c(pars$beta0_t), 
  iteration = rep(1:nrow(pars$beta0_t), times = ntree),
  tree      = factor(rep(1:ntree, each = nrow(pars$beta0_t))),
  plot      = factor(rep(plot_for_tree_idx, each = nrow(pars$beta0_t)))
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
      tree  = factor(1:ntree),
      plot  = factor(plot_for_tree_idx)
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
  parameter = 1:dim(beta_post)[2],
  knot      = 1:dim(beta_post)[3]
)
as.data.frame.table(beta_post, responseName = "beta") %>%
  mutate(
    iteration = factor(iteration), 
    parameter = factor(parameter), 
    knot      = factor(knot)
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
  sigma_true = rep(sigma_beta, each = nrow(pars$sigma_beta)),
  par        = factor(rep(1:length(sigma_beta), each = nrow(pars$sigma_beta)))
) %>%
  ggplot(aes(x = par, y = sigma_beta)) +
  geom_violin() +
  geom_point(alpha = 0.1, position = "jitter") +
  geom_point(aes(x = par, y = sigma_true), color = "red")  +
  ggtitle("spline penalty variance estimate")



## Fitted vs. estimated functional responses

dat_X <- data.frame(
  X           = c(X1, X2, X3),
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
  effects[k, , 1] <- X1_bs %*% pars$beta[k, 1, ]
  effects[k, , 2] <- X2_bs %*% pars$beta[k, 2, ]
  effects[k, , 3] <- X3_bs %*% pars$beta[k, 3, ]
}

dat_effects <- as.data.frame.table(effects, responseName = "effect")
dat_effects$iteration <- as.numeric(dat_effects$iteration)
dat_effects$observation <- as.numeric(dat_effects$observation)
dat_effects <- dat_effects %>%
  left_join(dat_X)

dat_truth <- data.frame(
  effect      = c(
    X1_bs %*% beta[1, ],
    X2_bs %*% beta[2, ],
    X3_bs %*% beta[3, ]
  ),
  X           = c(X1, X2, X3),
  observation = 1:length(dat$y),
  parameter   = factor(rep(paste("var", 1:3, sep="-"), each = length(dat$y)))
  
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
  facet_wrap(~ parameter, ncol = 1) +
  ggtitle("estimate in grey, simulated trend in red")



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
