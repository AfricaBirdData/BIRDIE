library(BIRDIE)

rm(list = ls())

config <- BIRDIE::configPreambJAGS(2017, server = FALSE)


# Diagnose convergence (Rhat) ---------------------------------------------

rhat_df <- data.frame()

for(i in 1:length(config$species)){

    sp_code <- config$species[i]

    print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

    # Skip species if less than 5 suitable sites were detected during fitting model fitting
    error_file <- setSpOutFilePath("Less_5_sites", config, sp_code, ".txt")

    if(file.exists(error_file)){
        message(paste0("Less_5_sites_", sp_code, "_", config$years_ch, ".txt"))
        next
    }

    # Else proceed with diagnostics
    fit <- readRDS(setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))

    fit_stats <- BIRDIE:::processJAGSoutput(fit, DIC = FALSE, params.omit = NULL)

    rhats <- lapply(fit_stats$Rhat, function(x) sum(abs(x - 1) > 0.1)) %>%
        as.data.frame()

    rhats$sp <- sp_code
    rhats$nobs <- length(fit_stats$mean$mu_t)
    rhats$years <- config$years_ch

    rhat_df <- rbind(rhat_df, rhats)

}

write.csv(rhat_df,
          file.path(config$out_dir, paste0("abu_diag_", config$years_ch, ".csv")),
          row.names = FALSE)


# Explore a particular fit ------------------------------------------------

library(dplyr)
library(ggplot2)
library(bayesplot)

# Load rhat dataframe if not in the workspace already
rhat_df <- read.csv(file.path(config$out_dir, paste0("abu_diag_", config$years_ch, ".csv")))

# Set species
sp_code <- 4

# Load fit
fit <- readRDS(setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))

# Load counts
counts <- utils::read.csv(setSpOutFilePath("abu_model_data", config, sp_code, ".csv"))

# Find Rhats greater than 1.1 or smaller than 0.9
source("beta/diagnoseJAGSrhat.R")

diagnoseJAGSrhat(fit_stats, param = "beta")
diagnoseJAGSrhat(fit_stats, param = "G")
diagnoseJAGSrhat(fit_stats, param = "phi")
diagnoseJAGSrhat(fit_stats, param = "sig.alpha")

# Plot mixing for those specific sites/years
mcmc_trace(fit, pars = vars(contains("beta[46,")))
mcmc_trace(fit, pars = vars(contains("sig.alpha[28]")))
mcmc_trace(fit, pars = vars(contains("phi[28]")))
mcmc_trace(fit, pars = vars(matches("phi\\[52\\]")))
mcmc_trace(fit, pars = vars(matches("mu.beta\\[2\\,.*\\]")))
mcmc_trace(fit, pars = vars(matches("G\\[4.*\\,1\\]")))

mcmc_intervals(fit, pars = vars(contains("sig.eps[14]")))
mcmc_intervals(fit, pars = vars(contains("phi")))

mcmc_intervals(fit, pars = vars(matches("beta\\[.*\\,1\\]")))
mcmc_intervals(fit, pars = vars(matches("beta\\[46\\,.*\\]")))
mcmc_intervals(fit, pars = vars(matches("mu_beta\\[2\\,.*\\]")))

mcmc_intervals(fit, pars = vars(matches("mu_beta")))

mcmc_intervals(fit, pars = vars(matches("sig.zeta")))

mcmc_intervals(fit, pars = vars(matches("mu_t\\[46\\,.*\\]")))

# Create plots
pers_theme <- ggplot2::theme_bw()
p <- BIRDIE::plotSsm2ss(fit = fit_stats, ssm_counts = counts, linear = TRUE,
                        plot_options = list(pers_theme = pers_theme,
                                            colors = c("#71BD5E", "#B590C7")))

p_sel <- p[as.character(unique(counts$LocationCode)[c(46, 59)])]

plot(p_sel[[1]]$plot)
print(p_sel[[1]]$data[[1]], n = Inf)

plot(p_sel[[2]]$plot)
print(p_sel[[2]]$data[[1]], n = Inf)

p_sel <- p[[1]]
plot(p_sel$plot)

p_sel[[1]]$data[[1]] %>%
    print(n = Inf)

dd <- lapply(p, "[[", 2) %>%
    lapply(.,"[[", 1) %>%
    bind_rows()

dd %>%
    filter(season == 1) %>%
    ggplot(aes(group = factor(site))) +
    geom_line(aes(x = year, y = mu_est, col = factor(site)), alpha = 0.6) +
    scale_color_viridis_d(option = "H", guide = "none")

dd %>%
    filter(mu_est > 2000)

ss <- dd %>%
    filter(mu_est > 2000) %>%
    pull(site) %>%
    unique()

dd %>%
    filter(site == ss) %>%
    print(n = Inf)

which(unique(counts$LocationCode) == ss)
