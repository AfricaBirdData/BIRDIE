# Simulations of continuous-time SSM model

library(dplyr)
library(ggplot2)
library(rstan)

options(mc.cores = parallel::detectCores() - 2)
rstan_options(auto_write = TRUE)

rm(list = ls())

# Number of steps
N <- 100

# Initial state
x0 <- 10

# Initial change due to summer (this should be estimated)
beta0 <- 0

# Initial change due to winter (this should be estimated)
lambda0 <- 0

# Vector of summer/winter/other indicators
season <- rep(c("s", "w", "o"), length.out = N)
summer <- ifelse(season == "s", 1, 0)
winter <- ifelse(season == "w", 1, 0)

# Vector of time increments
dt <- rnbinom(N-1, size = 5, mu = 2)
hist(dt)

# Define process standard deviation parameters
sig_alpha <- 1
sig_beta <- 3
sig_lambda <- 0.5

# Define observation error standard deviation parameters
sig_s <- 1.5
sig_w <- 0.5
sig_o <- 1


# Simulate process --------------------------------------------------------

# Define state vector
X <- matrix(NA, nrow = 3, ncol = N)

# Define initial state vector
X[,1] <- matrix(c(x0, beta0, lambda0))

# Iterate through steps
for(t in seq_len(N-1)){

    # transition matrix
    theta <- rbind(c(1, summer[t], winter[t]),
                   c(0, 1, 0),
                   c(0, 0, 1))

    # update process
    X[,t+1] = theta %*% X[,t,drop=FALSE] + matrix(rnorm(3, 0, c(sig_alpha, sig_beta, sig_lambda)*sqrt(dt[t])))
}

plot(X[1,])
plot(X[1, 1:20])

op <- par()
par(mfrow = c(3, 1))
plot(X[1,], type = "l")
plot(X[2,], type = "l")
plot(X[3,], type = "l")
par(op)


# Simulate multiple trajectories ------------------------------------------

# # Number of trajectories
# J <- 100
#
# # Define state vectors
# X <- array(NA, dim = c(3, N, J))
#
# # Define initial state vector
# for(j in seq_len(J)){
#     X[,1,j] <- matrix(c(x0, beta0, lambda0))
# }
#
# # Iterate through sims and through steps
# for(j in seq_len(J)){
#     for(t in seq_len(N-1)){
#
#         # transition matrix
#         theta <- rbind(c(1, summer[t], winter[t]),
#                        c(0, 1, 0),
#                        c(0, 0, 1))
#
#         # update process
#         X[,t+1,j] = theta %*% X[,t,j,drop=FALSE] + matrix(rnorm(3, 0, c(sig_alpha, sig_beta, sig_lambda)*sqrt(dt[t])))
#     }
# }
#
# plot(X[1,,1], type = "l", ylim = c(-150, 150))
# for(j in seq_len(J-1)){
#     lines(X[1,,j+1], col = j+1)
# }



# Fit model to recover parameters -----------------------------------------


# Bundle data
data.bundle <- list(  # number of observations
    N = ncol(X),
    # number of observed states
    D = 1,
    # number of states
    M = nrow(X),
    # observed data
    y = t(X[1,,drop = F]),
    # season indicators
    summer = summer,
    winter = winter,
    # observation matrix
    H = matrix(c(1, 0, 0), nrow = 1),
    # time increments
    dt = c(0.1, dt),
    # Priors for rate parameters
    psig_alpha = 3,
    psig_beta = 3,
    psig_lambda = 3,
    # Priors for observation error
    psig_s = 2,
    psig_w = 2,
    psig_o = 2
)

# Define parameters to plot
paramToPlot <- c("sig_alpha", "sig_beta", "sig_lambda",
                 "sig_s", "sig_w", "sig_o")

param <- c(paramToPlot)

# Compile model
stan_mod <- stan_model(file = "analysis/models/cwac_kalman_smooth.stan")

# Define initial values
init = function() list(
  mu0 = rnorm(1, 0, 1),
  beta0  = rnorm(1, 0, 1),
  lambda0  = rnorm(1, 0, 1),
  sig_alpha0 = rexp(1, 1),
  sig_beta0 = rexp(1, 1),
  sig_lambda0 = rexp(1, 1),

  sig_alpha = rexp(1, 1),
  sig_beta = rexp(1, 1),
  sig_lambda = rexp(1, 1),

  sig_s = rexp(1, 1),
  sig_w = rexp(1, 1),
  sig_o  = rexp(1, 1)
)

# gc()
start_time <- Sys.time()
fit <- sampling(stan_mod,
                data = data.bundle,
                pars = param, init = init,
                iter = 2000, chains = 4, warmup = 1000, thin = 1,
                cores = getOption("mc.cores", 4),
                control = list(adapt_delta = 0.85),
                refresh = max(2000/50,1),
                save_warmup = FALSE)
end_time <- Sys.time()
end_time - start_time

stan_trace(fit, pars = paramToPlot)
print(fit, pars = paramToPlot, probs = c(0.025, 0.5, 0.975),
      digits_summary = 3, use_cache = F)

