library(BIRDIE)
library(tidyverse)

rm(list = ls())

counts <- barberspan

ssmcounts <- prepCtSsmData(counts, species = 4)

fit <- fitCwacCtSsm(ssmcounts, mod_file = "analysis/models/cwac_ct_ssm_dyn.jags",
                    param = c("beta", "lambda", "sig.B", "sig.zeta", "sig.eps",
                              "sig.alpha", "sig.e", "sig.o", "mu"),
                    jags_control = list(ncores = 3))

summary(fit)
fit

out <- plotCtSsm(fit, ssmcounts)

plot(out$plot)

out$data


# Kalman smoother ---------------------------------------------------------

library(rstan)

options(mc.cores = parallel::detectCores() - 2)
rstan_options(auto_write = TRUE)

counts <- barberspan

ssmcounts <- prepCtSsmData(counts, species = 83)

# Bundle data
data.bundle <- list(  # number of observations
    N = nrow(ssmcounts),
    # number of observed states
    D = 1,
    # number of states
    M = 3,
    # observed data
    y = matrix(log(ssmcounts$count+1), ncol = 1),
    # season indicators
    summer = ssmcounts$summer,
    winter = ssmcounts$winter,
    # observation matrix
    H = matrix(c(1, 0, 0), nrow = 1),
    # time increments
    dt = ssmcounts$dt[!is.na(ssmcounts$dt)],
    # Priors for initial values
    pmu_mu0 = 5,
    psig_mu0 = 3,
    psig_beta0 = 3,
    psig_lambda0 = 3,
    # Priors for rate parameters
    psig_alpha = 2,
    psig_beta = 2,
    psig_lambda = 2,
    # Priors for observation error
    psig_s = 5,
    psig_w = 5,
    psig_o = 5
)

# Define parameters to plot
paramToPlot <- c("sig_alpha", "sig_beta", "sig_lambda",
                 "sig_s", "sig_w", "sig_o", "mu0", "beta0", "lambda0")

param <- c(paramToPlot, "pred")

# Compile model
stan_mod <- stan_model(file = "analysis/models/cwac_kalman_smooth.stan")

# Define initial values
init = function() list(
    mu0 = rexp(1, 0.1),
    beta0  = rnorm(1, 3, 3),
    lambda0  = rnorm(1, -3, 3),
    sig_alpha = rexp(1, 1),
    sig_beta = rexp(1, 1),
    sig_lambda = rexp(1, 1),
    sig_s = rexp(1, 1),
    sig_w = rexp(1, 1),
    sig_o  = rexp(1, 1)
)

fit <- sampling(stan_mod,
                data = data.bundle,
                pars = param, init = init,
                iter = 2000, chains = 4, warmup = 1000, thin = 1,
                cores = getOption("mc.cores", 4),
                control = list(adapt_delta = 0.95),
                refresh = max(2000/50,1),
                save_warmup = FALSE)

stan_trace(fit, pars = paramToPlot)
print(fit, pars = paramToPlot, probs = c(0.025, 0.5, 0.975),
      digits_summary = 3, use_cache = F)


musamples <- rstan::extract(fit, pars = "pred")$pred[,,1]

muhat <- data.frame(est = apply(musamples, 2, mean),
                    lb = apply(musamples, 2, quantile, 0.025),
                    ub = apply(musamples, 2, quantile, 0.975),
                    obs = data.bundle$y,
                    season = ssmcounts$season,
                    t = cumsum(c(0, data.bundle$dt)))

ggplot(muhat) +
    geom_point(aes(x = t, y = obs), col = "red") +
    geom_point(aes(x = t, y = est)) +
    geom_line(aes(x = t, y = lb), linetype = 2) +
    geom_line(aes(x = t, y = ub), linetype = 2) +
    facet_wrap("season", nrow = 3)

ggplot(muhat) +
    geom_point(aes(x = t, y = obs, col = season)) +
    geom_point(aes(x = t, y = est)) +
    geom_line(aes(x = t, y = est)) +
    geom_line(aes(x = t, y = lb), linetype = 2) +
    geom_line(aes(x = t, y = ub), linetype = 2)
