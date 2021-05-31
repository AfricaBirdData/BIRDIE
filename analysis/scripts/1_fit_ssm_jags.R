# 17-05-2021

rm(list = ls())

library(tidyverse)
library(CWAC)
library(BIRDIE)
library(jagsUI)


# Load data ---------------------------------------------------------------

# Barberspan data is already in the package

counts <- barberspan


# Prepare data to fit an SSM ----------------------------------------------

ssmcounts <- prepSsmData(counts, species = NULL)


# Fit 2-season fixed trend model ------------------------------------------

fit_fxd <- fitCwacSsm2ss(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))


# Fit 2-season dynamic trend model ----------------------------------------

fit_dyn <- fitCwacSsm2ss(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_dyn.jags",
                         param = c("beta", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))


# Explore models ----------------------------------------------------------

fit_fxd

fit_dyn

plotSsm2ss(fit = fit_fxd, ssm_counts = ssmcounts)

plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, dyn = TRUE)


# Fit model for resident species ------------------------------------------

fit_fxd <- fitCwacSsm2ss(ssmcounts, mod_file = "analysis/models/cwac_ssm_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t"))

fit_fxd

# Plot
plotSsm(fit = fit_fxd, ssm_counts = ssmcounts)

