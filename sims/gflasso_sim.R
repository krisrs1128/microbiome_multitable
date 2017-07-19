#! /usr/bin/env Rscript

# File description -------------------------------------------------------------

## ---- sparse-sim-packages ----
# List of packages for session
.packages = c("gflasso",
              "plyr",
              "dplyr",
              "reshape2",
              "glmnet",
              "Rcpp",
              "ggplot2")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.packages[!.inst], repos='http://cran.rstudio.com/')
}

# Load packages into session
lapply(.packages, require, character.only=TRUE)
set.seed(04032016)

cat("\014")  # Clear console

rm(list=ls()) # Delete all existing variables
graphics.off() # Close all open plots
theme_set(theme_bw())
sim_methods_dir <- "~/Documents/programming/simFuns/simMethods/"

## ---- simdata-params ----
K <- 200 # number of columns of Y
J <- 30 # number of columns of X
J0 <- 5 # number of columns of X actually contribute
n <- 40 # number of samples

## ---- gflasso-simulate-data ----
beta <- rnorm(J, .5, 1)
beta[sample(J)[1:(J / 2)]] <- rnorm(J / 2, -.5, 1)
B <- beta %*% t(rep(1, K)) + matrix(rnorm(J * K, 0, .1), J, K)
B[-sample(J, J0), 1 : (K / 2)] <- 0
B[-sample(J, J0), (K / 2 + 1) : K] <- 0
X <- matrix(rnorm(n * J), n, J)
Y <- X %*% B + matrix(rnorm(n * K, 0, .5), n, K)

## ---- vis-B ----
m_b <- melt(B, varnames = c("feature", "task"), value.name = "beta")

## ---- vis-reg-data ----
m_y <- melt(Y, varnames = c("sample", "task"), value.name = "y")
m_x <- melt(X, varnames = c("sample", "feature"), value.name = "x")
reg_data <- m_y %>%
  left_join(m_x) %>%
  left_join(m_b)
cur_data <- reg_data %>%
  filter(task %in% 95:105,
         feature %in% 1:10)

## ---- lasso-baseline ----
b_lasso <- matrix(0, J, K)
for (k in seq_len(K)) {
  glmnet_fit <- cv.glmnet(x = X, y = Y[, k], intercept = F)
  b_lasso[, k] <- coef(glmnet_fit)[-1]
}

## ---- gflasso ----
R <- matrix(1, K, K)
R[(K / 2 + 1) : K, 1:(K / 2)] <- 0
R[1:(K / 2), (K / 2 + 1) : K] <- 0
n_iter <- 200

lambda <- .15
gamma <- .05
eps <- .1
gflasso_res <- gflasso(Y, X, R,
                       list(delta_conv = 1, lambda = lambda, gamma = gamma,
                            iter_max = n_iter, verbose = T, eps = eps))

save(X, Y, B, b_lasso, gflasso_res, file = "gflasso_sim.RData")
