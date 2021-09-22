# 6-07-2021

rm(list = ls())

library(BIRDIE)


# Run SSM for one species -------------------------------------------------

# Load data (Barberspan example)
counts <- barberspan

# Select a species (ADU code)
sp <- 6

# Prepare data to fit an SSM
ssmcounts <- BIRDIE::prepSsmData(counts, species = sp)

# Fit 2-season dynamic trend model
fit_dyn <- BIRDIE::fitCwacSsm(ssmcounts, mod_file = BIRDIE::writeJagsModelFile(),
                              param = c("beta", "lambda", "sig.zeta", "sig.w",
                                        "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))

# Plot
pers_theme <- ggplot2::theme_bw()
p <- BIRDIE::plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, linear = FALSE,
                        plot_options = list(pers_theme = pers_theme,
                                            colors = c("#71BD5E", "#B590C7")))

plot(p$plot)


# Run SSM for all species -------------------------------------------------

# Load data (Barberspan example)
counts <- barberspan

# TAKES LONG
loopSsmAllSpp(barberspan, data_outdir = "analysis/out_nosync/",
              plot_outdir = "analysis/out_nosync/",
              param = c("beta", "lambda", "sig.zeta", "sig.w", "sig.eps",
                        "sig.alpha", "sig.e", "mu_t", "mu_wt"),
              jags_control = list(ncores = 3))
