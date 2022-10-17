library(BIRDIE)

rm(list = ls())

config <- BIRDIE::configPreambJAGS(2017, server = FALSE)


# Diagnose convergence (Rhat) ---------------------------------------------

rhat_df <- data.frame()

for(s in 1:length(config$species)){

    sp_code <- config$species[s]

    print(paste0("Working on species ", sp_code, " (", s, " of ", length(config$species), ")"))

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
diagnoseJAGSrhat(fit_stats, param = "sig.eps")

# Plot mixing for those specific sites/years
mcmc_trace(fit, pars = vars(contains("beta[46,")))
mcmc_trace(fit, pars = vars(contains("sig.e[31]")))
mcmc_trace(fit, pars = vars(contains("phi[28]")))
mcmc_trace(fit, pars = vars(matches("phi\\[41\\]")))
mcmc_trace(fit, pars = vars(matches("mu.beta\\[2\\,.*\\]")))
mcmc_trace(fit, pars = vars(matches("G\\[4.*\\,1\\]")))
mcmc_trace(fit, pars = vars(matches(".*\\[33\\]")))
mcmc_trace(fit, pars = vars(matches("lambda\\[33")))

mcmc_intervals(fit, pars = vars(contains("sig.zeta")))
mcmc_intervals(fit, pars = vars(contains("phi")))
mcmc_intervals(fit, pars = vars(matches("G\\[4.*\\,3\\]")))
mcmc_intervals(fit, pars = vars(matches("beta\\[.*\\,1\\]")))
mcmc_intervals(fit, pars = vars(matches("beta\\[34\\,.*\\]")))
mcmc_intervals(fit, pars = vars(matches("mu_beta\\[2\\,.*\\]")))

mcmc_intervals(fit, pars = vars(matches("mu_beta")))

mcmc_intervals(fit, pars = vars(matches("sig.zeta")))

mcmc_intervals(fit, pars = vars(matches("mu_t\\[46\\,.*\\]")))

# Create plots
pers_theme <- ggplot2::theme_bw()
p <- BIRDIE::plotSsm2ss(fit = fit_stats, ssm_counts = counts, linear = TRUE,
                        plot_options = list(pers_theme = pers_theme,
                                            colors = c("#71BD5E", "#B590C7")))

which(names(p) == 25452752)

p_sel <- p[[34]] # check 15 at least, also check 20 and 10 for blown-up uncertainty in final years, 34 for spikes at the beginning

plot(p_sel$plot)
print(p_sel$data[[1]], n = Inf)
print(p_sel$data[[2]], n = Inf)

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


names(p[[1]])

sub_data <- do.call("rbind", lapply(p, function(x){
    x$data[[1]][1:10,]
}))

sub_data %>%
    filter(year == 1996, season == 1) %>%
    filter(!is.na(count)) %>%
    print(n = Inf)

counts %>%
    filter(LocationCode == 25452752, Season == "S") %>%
    ggplot() +
    geom_line(aes(x = year, y = pdsi_mean))


counts %>%
    dplyr::group_by(LocationCode) %>%
    dplyr::mutate(site = dplyr::cur_group_id()) %>%
    dplyr::ungroup()
