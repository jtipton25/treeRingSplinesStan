# spatially-varying coefficients example
library(tidyverse)
library(splines)
library(treeRingSplines)
library(spam)

if (!dir.exists(here::here("images")))
  dir.create(here::here("images"))
if (!dir.exists(here::here("images", "spatially-varying")))
  dir.create(here::here("images", "spatially-varying"))

# number of plots
n_plot <- 10 # square root
p <- 3
df <- 20
n_per_tree <- 40
n_tree_per_plot <- 8

age <- rep(seq(0, 1, length.out = n_per_tree), n_tree_per_plot * n_plot^2)

climate_normal <- expand.grid(seq(0, 1, length = n_plot), seq(0, 1, length = n_plot))
# expand the climate normal with replicates
climate_normal <- sapply(1:ncol(climate_normal), function(i) {rep(climate_normal[, i], each = n_per_tree * n_tree_per_plot)}, simplify = "matrix")
colnames(climate_normal) <- c("temp", "precip")
climate_normal <- as.data.frame(climate_normal)

W <- cbind(1, build_spline_interactions(cbind(climate_normal[, 1], climate_normal[, 2]), df = df)[1, , ])
Q <- make_Q(df, 0.99)
beta0 <- c(rnorm(1), rmvnorm.prec(1, rep(0, df^2), 100 * Q))
beta <- matrix(0, df^2 + 1, p)
for (j in 1:p) {
  beta[, j] <- c(rnorm(1), rmvnorm.prec(1, rep(0, df^2), 100 * Q))
}
# simulate p covariates that vary by year
X <- matrix(rnorm(n_plot^2 * n_per_tree * p), n_plot^2 * n_per_tree, p)
age_idx  <- rep(1:n_per_tree, n_tree_per_plot * n_plot^2)
tree_idx <- rep(1:(n_tree_per_plot * n_plot^2), each = n_per_tree)
site_idx <- rep(1:(n_plot^2), each = n_per_tree * n_tree_per_plot)
X <- matrix(0, n_plot^2 * n_per_tree * n_tree_per_plot, p)
# super inelegant, but I think this works
for (i in 1:n_per_tree) {
  for (j in 1:(n_plot^2)) {
    for (k in 1:p) {
      X[age_idx == i & site_idx == j, k] <- rnorm(1)
    }
  }
}

# expand X with replicates
# X <- sapply(1:ncol(X), function(i) {rep(X[, i], each = n_tree_per_plot)})
# ~ bs(precip, df = df, intercept = FALSE) * bs(temp, df = df, intercept = FALSE))

## add in age effect
Q_age <- make_Q_1d(df, 0.99)
# beta_age <- drop(rmvnorm.prec(1, rep(0, df), Q_age))
beta_age <- -3 + rnorm(df, seq(1, 0, length.out = df), 0.1)
W_age <- build_spline(as.matrix(age), df)[1, , ]

## no plot random effects
sigma2 <- 0.1^2
sigma <- sqrt(sigma2)
y <- W %*% beta0 + W_age %*% beta_age + apply(X * (W %*% beta), 1, sum) + rnorm(n_plot^2 * n_per_tree * n_tree_per_plot, 0, sigma)


p0 <- data.frame(temp = climate_normal$temp, precip = climate_normal$precip, beta0 = W %*% beta0) %>%
  ggplot(aes(x = temp, y = precip, fill = beta0)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("beta0 in climate space")
p1 <- data.frame(temp = climate_normal$temp, precip = climate_normal$precip, beta = W %*% beta[, 1]) %>%
  ggplot(aes(x = temp, y = precip, fill = beta)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("beta1 in climate space")
p2 <- data.frame(temp = climate_normal$temp, precip = climate_normal$precip, beta = W %*% beta[, 2]) %>%
  ggplot(aes(x = temp, y = precip, fill = beta)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("beta2 in climate space")
# p3 <- data.frame(temp = climate_normal$temp, precip = climate_normal$precip, beta = W %*% beta[, 3]) %>%
#   ggplot(aes(x = temp, y = precip, fill = beta)) +
#   geom_raster() +
#   scale_fill_viridis_c() +
#   ggtitle("beta3 in climate space")
p0 + p1 + p2
# (p0 + p1) / (p2 + p3)

dat <- data.frame(y = y,
           age = age,
           tree = factor(rep(1:(n_tree_per_plot * n_plot^2), each = n_per_tree)),
           site = factor(rep(1:(n_plot^2), each = n_per_tree * n_tree_per_plot)))

dat %>%
  # filter(site == 1) %>%
  filter(site %in% 1:20) %>%
  ggplot(aes(x = age, y = y, group = tree, color = tree)) +
  geom_line(alpha = 0.1) +
  facet_wrap(~ site, ncol = 10) +
  theme(legend.position = "none")

# data.frame(temp = climate_normal$temp, precip = climate_normal$precip, y = y) %>%
#   ggplot(aes(x = temp, y = precip, fill = y)) +
#   geom_raster() +
#   scale_fill_viridis_c() +
#   ggtitle("observations in climate space")


params <- list(
  n_adapt   = 200,
  n_mcmc    = 200,
  n_message = 5,
  n_thin    = 1)

out <- mcmc_svc_age(
  y               = drop(y),
  X               = X,
  age             = age,
  climate_normal  = as.matrix(climate_normal),
  params          = params,
  df              = df)

str(out)


## trace plot for intercept, age effect, and regression parameters
layout(matrix(1:4, 2, 2))
matplot(out$beta0, type = 'l')
matplot(out$beta_age, type = 'l')
matplot(out$beta[, , 1], type = 'l')
matplot(out$beta[, , 2], type = 'l')

## plot intercept, age effect, and residual fitted vs. estimated
layout(matrix(1:4, 2, 2))
plot(beta0, apply(out$beta0, 2, mean))
abline(0, 1, col = "red")

plot(beta_age, apply(out$beta_age, 2, mean))
abline(0, 1, col = "red")

hist(out$sigma2)
abline(v = sigma2, col = "red", lwd = 2)

## plot penalty parameter trace plots
layout(matrix(1:4, 2, 2))
hist(out$tau20)
abline(v = 100, col = "red")

for (j in 1:3) {
  hist(out$tau2[, j])
  abline(v = 100, col = "red")
}

# layout(matrix(1:4, 2, 2))
# plot(beta_age, apply(out$beta_age, 2, mean))
# abline(0, 1, col = "red")
#
# hist(out$sigma2)
# abline(v = sigma2, col = "red", lwd = 2)



layout(matrix(1:4, 2, 2))
for (j in 1:p) {
  plot(beta[, j], apply(out$beta[, , j], 2, mean))
  abline(0, 1, col = "red")
}

# trace plots for penalty parameters
layout(matrix(1:4, 2, 2))
plot(out$tau20, type = 'l')
plot(out$tau2_age, type = 'l')
matplot(out$tau2, type = 'l')


# plot fitted vs. observed beta0
data.frame(
  observation = 1:(n_plot^2),
  temp        = climate_normal$temp,
  precip      = climate_normal$precip,
  observed    = c(W %*% beta0),
  fitted      = c(out$W %*% apply(out$beta0, 2, mean)),
  parameter   = rep("beta0", each = n_plot^2)
) %>%
  pivot_longer(cols = c("fitted", "observed"), names_to = "type") %>%
  ggplot(aes(x = temp, y = precip, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(parameter ~ type) +
  ggtitle("Intercept")

# plot fitted vs. observed beta
p_fitted <- data.frame(
  observation = 1:(n_plot^2 * n_per_tree * n_tree_per_plot),
  temp        = climate_normal$temp,
  precip      = climate_normal$precip,
  observed    = c(W %*% beta),
  fitted      = c(out$W %*% apply(out$beta, c(2, 3), mean)),
  parameter   = rep(paste0("beta[", 1:p, "]"), each = n_plot^2 * n_per_tree * n_tree_per_plot)
) %>%
  group_by(parameter) %>%
  mutate(observed = (observed - mean(observed)) / sd(observed),
         fitted   = (fitted - mean(fitted)) / sd(fitted)) %>%
  ungroup() %>%
  pivot_longer(cols = c("fitted", "observed"), names_to = "type") %>%
  ggplot(aes(x = temp, y = precip, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(parameter ~ type) +
  ggtitle("Varying slope estimates (centered and scaled)")

p_fitted


ggsave(
  file = here::here("images", "spatially-varying", "p_fitted.png"),
  width = 16,
  height = 9,
  p_fitted
)
## next steps -- put this into tree ring framework

