# 17-05-2021

rm(list = ls())

library(tidyverse)
library(CWAC)
library(BIRDIE)
library(jagsUI)

figuredir <- "analysis/output/barberspan_plots/"


# Load data ---------------------------------------------------------------

# Barberspan data is already in the package

counts <- barberspan

sort(table(counts$spp))


# Fixed trend all species -------------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = NULL)

# Fit 2-season fixed trend model
fit_fxd <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"),
                         jags_control = list(ncores = 3))

# Plot
p <- plotSsm2ss(fit = fit_fxd, ssm_counts = ssmcounts)

ggsave(paste0(figuredir, "ssm_fxd_all.png"), plot = p$comb)


# Fixed trend common species -------------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = 212)

# Fit 2-season fixed trend model
fit_fxd <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))

# Plot
p <- plotSsm2ss(fit = fit_fxd, ssm_counts = ssmcounts)

ggsave(paste0(figuredir, "ssm_fxd_common.png"), plot = p$comb)


# Fixed trend rare species -------------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = 104)

# Fit 2-season fixed trend model
fit_fxd <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))

# Plot
p <- plotSsm2ss(fit = fit_fxd, ssm_counts = ssmcounts)

ggsave(paste0(figuredir, "ssm_fxd_rare.png"), plot = p$comb)


# Dynamic trend all species -------------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = NULL)

# Fit 2-season dynamic trend model
fit_dyn <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_dyn.jags",
                         param = c("beta", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))

# Plot
p <- plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, dyn = TRUE)

ggsave(paste0(figuredir, "ssm_dyn_all.png"), plot = p$comb)


# Dynamic trend common species -------------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = 212)

# Fit 2-season dynamic trend model
fit_dyn <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_dyn.jags",
                         param = c("beta", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))

# Plot
p <- plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, dyn = TRUE)

ggsave(paste0(figuredir, "ssm_dyn_common.png"), plot = p$comb)


# Dynamic trend rare species -------------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = 104)

# Fit 2-season dynamic trend model
fit_dyn <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_dyn.jags",
                         param = c("beta", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))

# Plot
p <- plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, dyn = TRUE)

ggsave(paste0(figuredir, "ssm_dyn_rare.png"), plot = p$comb)


# Fixed trend all species one season ---------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = NULL)

# Fit 1-season fixed trend model
fit_1ss <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t"))

# Plot
p <- plotSsm(fit = fit_1ss, ssm_counts = ssmcounts)

ggsave(paste0(figuredir, "ssm_1ss_all.png"), plot = p$comb)


# Fixed trend one season common species -------------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = 212)

# Fit 1-season fixed trend model
fit_1ss <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t"))

# Plot
p <- plotSsm(fit = fit_1ss, ssm_counts = ssmcounts)

ggsave(paste0(figuredir, "ssm_1ss_common.png"), plot = p$comb)


# Fixed trend one season rare species -------------------------------------------------

# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = 104)

# Fit 1-season fixed trend model
fit_1ss <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_fxd.jags",
                         param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t"))

# Plot
p <- plotSsm(fit = fit_1ss, ssm_counts = ssmcounts)

ggsave(paste0(figuredir, "ssm_1ss_rare.png"), plot = p$comb)




# Extract plot data -------------------------------------------------------


# Prepare data to fit an SSM
ssmcounts <- prepSsmData(counts, species = NULL)

# Fit 2-season dynamic trend model
fit_dyn <- fitCwacSsm(ssmcounts, mod_file = "analysis/models/cwac_ssm_2ss_dyn.jags",
                      param = c("beta", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"),
                      jags_control = list(ncores = 3))

# Plot
p <- plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, dyn = TRUE)

# Extract plot data
plotdata <- lapply(p$ind, ggplot_build)
plotdata <- lapply(plotdata, "[[", "data")

# Export sample
write.csv(plotdata[[1]][[1]], "analysis/output/plotdata_test.csv")

ggsave(paste0(figuredir, "ssm_dyn_all.png"), plot = p$comb)
