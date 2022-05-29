#' Summarize predictions from state-space JAGS model
#'
#' @param sp_code SAFRING code of the species of interest
#' @inheritParams ppl_run_pipe_abu1
#'
#' @return
#' @export
#'
#' @examples
ppl_summarize_ssm <- function(sp_code, config, ...){

    print(paste("Summarizing state-space JAGS model at", Sys.time()))

    # Load counts
    counts <- utils::read.csv(file.path(config$out_dir, sp_code, paste0("abu_model_data_", sp_code, "_", config$years_ch,".csv")))


    fit <- readRDS(file.path(config$out_dir, sp_code, paste0("ssm_fit_", config$years_ch, "_", sp_code, ".rds")))

    # Plot
    pers_theme <- ggplot2::theme_bw()
    p <- BIRDIE::plotSsm2ss(fit = fit, ssm_counts = counts, linear = TRUE,
                            plot_options = list(pers_theme = pers_theme,
                                                colors = c("#71BD5E", "#B590C7")))

    # grid::grid.newpage()
    # grid::grid.draw(p$plot)

    out_all_sites <- vector("list", length(unique(counts$LocationCode)))

    for(i in seq_along(unique(counts$LocationCode))){

        site_sel <- names(p)[[i]]

        plotdir <- file.path(config$out_dir, sp_code, paste0("ssm_plot_", sp_code, "_", site_sel, "_", config$years_ch, ".png"))
        ggplot2::ggsave(plotdir, plot = p[[i]]$plot, width = 6.4, height = 5.2)

        # Extract plot data
        plotdata <- p[[i]]$data

        # Prepare state data
        stt_df <- plotdata[[1]]

        stt_dfw <- cbind(stt_df %>%
                             dplyr::filter(season == 1) %>%
                             dplyr::select(-season),
                         stt_df %>%
                             dplyr::filter(season == 2) %>%
                             dplyr::select(-c(year, season, site)))

        names(stt_dfw) <- c("year", "site", "summer.count", "summer.est", "summer.ci.lower", "summer.ci.upper",
                            "winter.count", "winter.est", "winter.ci.lower", "winter.ci.upper")

        # Prepare trend data
        trd_df <- plotdata[[2]]

        names(trd_df) <- c("slope.est", "slope.ci.lower", "slope.ci.upper",
                           "winter.prop.est", "winter.prop.ci.lower", "winter.prop.ci.upper",
                           "year")

        # Combine
        out_df <- cbind(stt_dfw, dplyr::select(trd_df, -year))

        # Add missing columns
        out_df <- out_df %>%
            dplyr::mutate(species = sp_code)#,
        # summer.count = round(exp(log.summer.count)),
        # winter.count = round(exp(log.winter.count)))

        # Order columns
        out_df <- out_df %>%
            dplyr::select(species, site, year, summer.count, winter.count,
                          summer.est, winter.est, winter.ci.lower,
                          winter.ci.upper,
                          summer.ci.lower, summer.ci.upper,
                          slope.est, slope.ci.lower, slope.ci.upper,
                          winter.prop.est, winter.prop.ci.lower,
                          winter.prop.ci.upper)

        out_all_sites[[i]] <- out_df

    }

    out_all_sites <- dplyr::bind_rows(out_all_sites)

    # Export sample
    datadir <- file.path(config$out_dir, sp_code, paste0("ssm_pred_", sp_code, "_", site, "_", config$years_ch, ".csv"))
    utils::write.csv(out_df, datadir, row.names = FALSE)

}
