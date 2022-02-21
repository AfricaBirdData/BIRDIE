#' Summarize predictions from state-space JAGS model
#'
#' @inheritParams ppl_run_pipe_abu
#'
#' @return
#' @export
#'
#' @examples
ppl_summarize_ssm <- function(sp_code, site, year, config, ...){

    print(paste("Summarizing state-space JAGS model at", Sys.time()))

    # Load counts
    counts <- readRDS(paste0(config$data_dir, site, "_93_20_visit_covts.rds"))

    # Filter years
    counts <- counts %>%
        dplyr::filter(Year %in% config$years)

    # Prepare data to fit an SSM
    ssmcounts <- prepSsmData(counts, sp_code, keep = Hmisc::Cs(CountCondition, prcp, tmmn, tmmx, watext, watrec))

    # Load fit
    if(length(sp_code) > 1){
        sp_code <- "group"
    }

    fit <- readRDS(file.path(config$data_outdir, sp_code, paste0("ssm_fit_", site, "_", config$years_ch, "_", sp_code, ".rds")))

    # Plot
    pers_theme <- ggplot2::theme_bw()
    p <- BIRDIE::plotSsm2ss(fit = fit, ssm_counts = ssmcounts, linear = TRUE,
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
        dplyr::mutate(species = sp_code,
                      summer.count = round(exp(log.summer.count)),
                      winter.count = round(exp(log.winter.count)))

    # Order columns
    out_df <- out_df %>%
        dplyr::mutate(site = site) %>%
        dplyr::select(site, species, year, summer.count, winter.count,
                      log.summer.count, log.winter.count,
                      summer.logest, winter.logest, winter.logest.ci.lower,
                      winter.logest.ci.upper,
                      summer.logest.ci.lower, summer.logest.ci.upper,
                      slope.est, slope.ci.lower, slope.ci.upper,
                      winter.prop.est, winter.prop.ci.lower,
                      winter.prop.ci.upper)

    # Export sample
    datadir <- paste0(config$data_outdir, sp_code, "/ssm_pred_", site, "_", config$years_ch, "_", sp_code, ".csv")
    utils::write.csv(out_df, datadir, row.names = FALSE)

    plotdir <- paste0(config$plot_outdir, sp_code, "/ssm_plot_", site, "_", config$years_ch, "_", sp_code, ".png")
    ggplot2::ggsave(plotdir, plot = p$plot, width = 6.4, height = 5.2)

}
