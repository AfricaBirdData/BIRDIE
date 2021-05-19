# 17-05-2021

rm(list = ls())

library(tidyverse)
library(CWAC)
library(jagsUI)


# Load data ---------------------------------------------------------------

counts <- readRDS("analysis/data/barberspan.rds")

names(counts)


# Prepare season info -----------------------------------------------------

unique(counts$Season)

# Some records are classified as "O" (other). We are only interested in summer and winter
# so we filter out "O"
counts <- counts %>%
  filter(Season != "O")

# We also create a dummy season variable
counts <- counts %>%
  mutate(season_id = if_else(Season == "W", 2, 1))


# Filter migrant options --------------------------------------------------

# Migrant options
unique(counts$Migrant)

# Remove those species with no migrant status
counts %>%
  filter(Migrant == "")

counts <- counts %>%
  filter(Migrant != "")

resdt <- counts %>%
  filter(Migrant == "y")

migrt <- counts %>%
  filter(Migrant == "n")


# Fit resident JAGS model -------------------------------------------------

# Prepare data for resident species:

# Sum counts for each year across species
resdt <- resdt %>%
  mutate(year = lubridate::year(startDate)) %>%
  group_by(year, startDate, season_id) %>%
  summarize(count = sum(count)) %>%
  ungroup()

# Fill in years with no counts. Note that there are a number of years with
# no counts that will be treated as missing data
resdt <- resdt %>%
  complete(year = min(year):max(year), season_id)

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
