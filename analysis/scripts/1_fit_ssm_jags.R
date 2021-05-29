# 17-05-2021

rm(list = ls())

library(tidyverse)
library(CWAC)
library(BIRDIE)
library(jagsUI)


# Load data ---------------------------------------------------------------

# Find location code
loc_code <- listCwacSites("North West") %>%
  filter(Name == "Barberspan") %>%
  pull(Loc_code)

counts <- getCwacSiteCounts(loc_code)

names(counts)

counts <- barberspan

# Prepare data to fit an SSM ----------------------------------------------

ssmcounts <- prepSsmData(counts)


# Fit 2-season fixed trend model ------------------------------------------

fit_fxd <- fitCwacSsm2ss(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))


# Fit 2-season dynamic trend model ----------------------------------------

fit_dyn <- fitCwacSsm2ss(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_dyn.jags",
                         param = c("beta", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))


# Explore models ----------------------------------------------------------

fit_fxd

fit_dyn

plotSsm(fit = fit_fxd, ssm_counts = ssmcounts)


# Fit resident JAGS model -------------------------------------------------

# Prepare data for resident species

mod.data <- list(y = log(resdt$count),
                    n = nrow(resdt),
                    x = resdt$season_id)


# Inits function
mod.inits <- function(){
  list (#b = c(0, runif(1, log(1),log(1000))),
        tau.w2 = runif(1, log(1),log(5)),
        tau.eps2 = runif(1, log(1),log(5)),
        tau.alpha = runif(1, log(1),log(5)))
}


# MCMC settings
ni <- 10000 # number of iterations
nb <- 5000  # burning iterations
nt <- 1     # chain thinning
nc <- 3     # number of chains
na <- NULL   # (default) adapting iterations

# Parameters to estimate
mod.param <- c("mu_t", "sig.w2", "sig.eps2", "sig.alpha")

# Start Gibbs sampling
mod.fit <- jags(data = mod.data, inits = mod.inits,
                parameters.to.save = mod.param, model.file = "analysis/models/cwac_ssm_resident.jags",
                n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = 1,
                modules = c('glm','lecuyer', 'dic'), factories = NULL, parallel = T, n.cores = nc, DIC = TRUE,
                verbose = TRUE)


# Explore model for residents -------------------------------------------------

# Summary
summary(mod.fit)
mod.fit

traceplot(mod.fit, parameters = "mu_t[1:4]")

# Plot time series
post_stt_resdt <- data.frame(est = mod.fit$mean$mu_t,
                             lb = mod.fit$q2.5$mu_t,
                             ub = mod.fit$q97.5$mu_t) %>%
  mutate(year = resdt$year,
         season = resdt$season_id,
         count = log(resdt$count)) %>%
  group_by(year, season) %>%
  mutate(seas_id = cur_group_id()) %>%
  ungroup()

# In a single plot
post_stt_resdt %>%
  pivot_longer(cols = c(est, lb, ub),
               names_to = "quantile") %>%
  ggplot() +
  geom_path(aes(x = seas_id, y = value, linetype = quantile)) +
  geom_point(aes(x = seas_id, y = count, col = factor(season))) +
  scale_linetype_manual(values = c(1, 2, 2))

# Separated by season
post_stt_resdt %>%
  pivot_longer(cols = c(est, lb, ub),
               names_to = "quantile") %>%
  ggplot() +
  geom_path(aes(x = year, y = value, linetype = quantile)) +
  geom_point(aes(x = year, y = count, col = factor(season))) +
  scale_linetype_manual(values = c(1, 2, 2)) +
  facet_wrap("season", nrow = 2)


# Fit migrants JAGS model -------------------------------------------------

# Prepare data for resident species

mod.data <- list(winter = log(migrt %>% filter(season_id == 2) %>% pull(count)),
                 summer = log(migrt %>% filter(season_id == 1) %>% pull(count)),
                 N = nrow(migrt)/2)


# Inits function
# mod.inits <- function(){
#   list (#b = c(0, runif(1, log(1),log(1000))),
#     tau.w2 = runif(1, log(1),log(5)),
#     tau.eps2 = runif(1, log(1),log(5)),
#     tau.alpha = runif(1, log(1),log(5)))
# }


# MCMC settings
ni <- 10000 # number of iterations
nb <- 5000  # burning iterations
nt <- 1     # chain thinning
nc <- 3     # number of chains
na <- NULL   # (default) adapting iterations

# Parameters to estimate
mod.param <- c("mu_t", "mu_wt", "lambda")

# Start Gibbs sampling
mod.fit <- jags(data = mod.data, #inits = mod.inits,
                parameters.to.save = mod.param, model.file = "analysis/models/cwac_ssm_migrant.jags",
                n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = 1,
                modules = c('glm','lecuyer', 'dic'), factories = NULL, parallel = T, n.cores = nc, DIC = TRUE,
                verbose = TRUE)


# Explore model for migrants -------------------------------------------------

# Summary
summary(mod.fit)
mod.fit

traceplot(mod.fit, parameters = "mu_t[1:4]")

# Plot time series
post_stt_resdt <- data.frame(est = c(mod.fit$mean$mu_t, mod.fit$mean$mu_wt),
                             lb = c(mod.fit$q2.5$mu_t, mod.fit$q2.5$mu_wt),
                             ub = c(mod.fit$q97.5$mu_t, mod.fit$q97.5$mu_wt),
                             year = rep(unique(migrt$year), 2),
                             season = rep(unique(migrt$season_id), each = mod.data$N),
                             count = c(mod.data$summer, mod.data$winter)) %>%
  group_by(year, season) %>%
  mutate(seas_id = cur_group_id()) %>%
  ungroup()


# Plot separated by season
post_stt_resdt %>%
  pivot_longer(cols = c(est, lb, ub),
               names_to = "quantile") %>%
  ggplot() +
  geom_path(aes(x = year, y = value, linetype = quantile)) +
  geom_point(aes(x = year, y = count, col = factor(season))) +
  scale_linetype_manual(values = c(1, 2, 2)) +
  facet_wrap("season", nrow = 2)


# Fit migrants JAGS model to resident species ---------------------------

# Prepare data for resident species

mod.data <- list(winter = log(resdt %>% filter(season_id == 2) %>% pull(count)),
                 summer = log(resdt %>% filter(season_id == 1) %>% pull(count)),
                 N = nrow(resdt)/2)


# Inits function
# mod.inits <- function(){
#   list (#b = c(0, runif(1, log(1),log(1000))),
#     tau.w2 = runif(1, log(1),log(5)),
#     tau.eps2 = runif(1, log(1),log(5)),
#     tau.alpha = runif(1, log(1),log(5)))
# }


# MCMC settings
ni <- 10000 # number of iterations
nb <- 5000  # burning iterations
nt <- 1     # chain thinning
nc <- 3     # number of chains
na <- NULL   # (default) adapting iterations

# Parameters to estimate
mod.param <- c("mu_t", "mu_wt", "lambda", "sig.zeta", "sig.w2", "sig.alpha", "sig.e", "tau.eps")

# Start Gibbs sampling
mod.fit2 <- jags(data = mod.data, #inits = mod.inits,
                parameters.to.save = mod.param, model.file = "analysis/models/cwac_ssm_migrant.jags",
                n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = 1,
                modules = c('glm','lecuyer', 'dic'), factories = NULL, parallel = T, n.cores = nc, DIC = TRUE,
                verbose = TRUE)


# Explore model for migrants -------------------------------------------------

# Summary
summary(mod.fit2)
mod.fit2

traceplot(mod.fit2, parameters = "mu_t[1:4]")

# Plot time series
post_stt_resdt <- data.frame(est = c(mod.fit2$mean$mu_t, mod.fit2$mean$mu_wt),
                             lb = c(mod.fit2$q2.5$mu_t, mod.fit2$q2.5$mu_wt),
                             ub = c(mod.fit2$q97.5$mu_t, mod.fit2$q97.5$mu_wt),
                             year = rep(unique(resdt$year), 2),
                             season = rep(unique(resdt$season_id), each = mod.data$N),
                             count = c(mod.data$summer, mod.data$winter)) %>%
  group_by(year, season) %>%
  mutate(seas_id = cur_group_id()) %>%
  ungroup()


# Plot separated by season
post_stt_resdt %>%
  pivot_longer(cols = c(est, lb, ub),
               names_to = "quantile") %>%
  ggplot() +
  geom_path(aes(x = year, y = value, linetype = quantile)) +
  geom_point(aes(x = year, y = count, col = factor(season))) +
  scale_linetype_manual(values = c(1, 2, 2)) +
  facet_wrap("season", nrow = 2)
