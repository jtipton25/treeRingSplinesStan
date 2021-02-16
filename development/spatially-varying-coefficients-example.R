# spatially-varying coefficients example
library(tidyverse)
library(splines)
library(treeRingSplines)
library(spam)

if (!dir.exists(here::here("images")))
  dir.create(here::here("images"))
if (!dir.exists(here::here("images", "spatially-varying")))
  dir.create(here::here("images", "spatially-varying"))

N <- 80
p <- 6
df <- 20
climate_normal <- expand.grid(seq(0, 1, length = N), seq(0, 1, length = N))
colnames(climate_normal) <- c("temp", "precip")
W <- cbind(1, build_spline_interactions(cbind(climate_normal[, 1], climate_normal[, 2]), df = df)[1, , ])
Q <- make_Q(df, 0.99)
beta0 <- c(rnorm(1), rmvnorm.prec(1, rep(0, df^2), Q))
beta <- matrix(0, df^2 + 1, p)
for (j in 1:p) {
  beta[, j] <- c(rnorm(1), rmvnorm.prec(1, rep(0, df^2), Q))
}
# simulate p covariates
X <- matrix(rnorm(N^2 * p), N^2, p)
# ~ bs(precip, df = df, intercept = FALSE) * bs(temp, df = df, intercept = FALSE))

sigma2 <- 0.1^2
sigma <- sqrt(sigma2)
y <- W %*% beta0 + apply(X * (W %*% beta), 1, sum) + rnorm(N, 0, sigma)


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

data.frame(temp = climate_normal$temp, precip = climate_normal$precip, y = y) %>%
  ggplot(aes(x = temp, y = precip, fill = y)) +
  geom_raster() +
  scale_fill_viridis_c() +
  ggtitle("observations in climate space")


params <- list(
  n_adapt   = 200,
  n_mcmc    = 200,
  n_message = 5,
  n_thin    = 1)

out <- mcmc_svc(
  y               = drop(y),
  X               = X,
  climate_normal  = as.matrix(climate_normal),
  params          = params,
  df              = df)

str(out)

layout(matrix(1:4, 2, 2))
plot(beta0, apply(out$beta0, 2, mean))
abline(0, 1, col = "red")

hist(out$sigma2)
abline(v = sigma2)


layout(matrix(1:4, 2, 2))
for (j in 1:p) {
  plot(beta[, j], apply(out$beta[, , j], 2, mean))
  abline(0, 1, col = "red")
}



# plot fitted vs. observed beta0
data.frame(
  observation = 1:(N^2),
  temp        = climate_normal$temp,
  precip      = climate_normal$precip,
  observed    = c(W %*% beta0),
  fitted      = c(out$W %*% apply(out$beta0, 2, mean)),
  parameter   = rep("beta0", each = N^2)
) %>%
  pivot_longer(cols = c("fitted", "observed"), names_to = "type") %>%
  ggplot(aes(x = temp, y = precip, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(parameter ~ type) +
  ggtitle("Intercept")

# plot fitted vs. observed beta
p_fitted <- data.frame(
  observation = 1:(N^2),
  temp        = climate_normal$temp,
  precip      = climate_normal$precip,
  observed    = c(W %*% beta),
  fitted      = c(out$W %*% apply(out$beta, c(2, 3), mean)),
  parameter   = rep(paste0("beta[", 1:p, "]"), each = N^2)
) %>%
  pivot_longer(cols = c("fitted", "observed"), names_to = "type") %>%
  ggplot(aes(x = temp, y = precip, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(parameter ~ type) +
  ggtitle("Varying slope estimates")


ggsave(
  file = here::here("images", "spatially-varying", "p_fitted.png"),
  width = 16,
  height = 9,
  p_fitted
)
## next steps -- put this into tree ring framework

