library(BIRDIE)

rm(list = ls())


# Script settings ---------------------------------------------------------

mod_file <- "/drv_birdie/birdie_ftp/Working/git/BIRDIE/analysis/models/cwac_ssm_2ss_dyn.jags"
data_outdir <- "/drv_birdie/birdie_ftp/"
plot_outdir <- "/drv_birdie/birdie_ftp/"


# Prepare count data ------------------------------------------------------

counts <- barberspan

# Identify site
site_id <- unique(counts$LocationCode)

if(length(site_id) != 1){
    stop("Either location code is missing or there are multiple sites.")
}

# Get species list
spp <- unique(counts$SppRef)


# Fit models and save results ---------------------------------------------

for(i in seq_along(spp)){

    sp <- spp[i]

    # Prepare data to fit an SSM
    ssmcounts <- BIRDIE::prepSsmData(counts, spp_sel = sp)

    # Fit 2-season dynamic trend model
    fit_dyn <- BIRDIE::fitCwacSsm(ssmcounts, mod_file = mod_file,
                                  param = c("beta", "lambda", "sig.zeta", "sig.w",
                                            "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"),
                                  jags_control = list(ncores = 3))

    # Plot
    pers_theme <- ggplot2::theme_bw()
    p <- BIRDIE::plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, linear = TRUE,
                            plot_options = list(pers_theme = pers_theme,
                                                colors = c("#71BD5E", "#B590C7")))

    # grid::grid.newpage()
    # grid::grid.draw(p$plot)

    # Extract plot data
    plotdata <- p$data

    # Prepare state data
    stt_df <- plotdata[[1]]

    stt_dfw <- cbind(stt_df %>%
                         dplyr::filter(season == 1) %>%
                         dplyr::select(-season),
                     stt_df %>%
                         dplyr::filter(season == 2) %>%
                         dplyr::select(-c(year, season)))

    names(stt_dfw) <- c("summer.logest", "summer.logest.ci.lower", "summer.logest.ci.upper",
                        "year", "log.summer.count",
                        "winter.logest", "winter.logest.ci.lower", "winter.logest.ci.upper",
                        "log.winter.count")

    # Prepare trend data
    trd_df <- plotdata[[2]]

    names(trd_df) <- c("slope.est", "slope.ci.lower", "slope.ci.upper",
                       "winter.prop.est", "winter.prop.ci.lower", "winter.prop.ci.upper",
                       "year")

    # Combine
    out_df <- cbind(stt_dfw, dplyr::select(trd_df, -year))

    # Add missing columns
    out_df <- out_df %>%
        dplyr::mutate(species = sp,
                      summer.count = round(exp(log.summer.count)),
                      winter.count = round(exp(log.winter.count)))

    # Order columns
    out_df <- out_df %>%
        dplyr::select(species, year, summer.count, winter.count,
                      log.summer.count, log.winter.count,
                      summer.logest, winter.logest, winter.logest.ci.lower,
                      winter.logest.ci.upper,
                      summer.logest.ci.lower, summer.logest.ci.upper,
                      slope.est, slope.ci.lower, slope.ci.upper,
                      winter.prop.est, winter.prop.ci.lower,
                      winter.prop.ci.upper)

    # Export sample
    datadir <- paste0(data_outdir, sp, "/ssm_dat_", site_id, "_", sp, ".csv")
    utils::write.csv(out_df, datadir, row.names = FALSE)

    plotdir <- paste0(plot_outdir, sp, "/ssm_plot_", site_id, "_", sp, ".png")
    ggplot2::ggsave(plotdir, plot = p$plot, width = 6.4, height = 5.2)

}
