#' Plot state series of a CWAC state-space model with two seasons
#'
#' @param fit A JAGS state-space model fitted to CWAC data
#' @param ssm_counts A data frame with the count data use to fit the
#' state-space model
#' @param linear If TRUE (default) abundance estimates and data are
#' transformed back to its orginal scale.
#' @param plot_options A list with two elements: colors - the colours of the
#' points that will appear in the plot (two values), and pers_theme  - A
#' personalized ggplot theme.
#'
#' @return A list with two elements: i) plot: a plot with summer and winter fitted states, as well as the long-term trend,
#' ii) data: the data used to create the individual plots. This is useful for extracting the data used by ggplot to render the plots
#' (e.g. for exporting to the dashboard)
#' @importFrom ggplot2 aes
#' @export
#'
#' @examples
#' counts <- barberspan
#' ssmcounts <- prepSsmData(counts, species = NULL)
#' fit <- fitCwacSsm(ssmcounts, param = c("beta", "lambda", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))
#' plotSsm2ss(fit = fit, ssm_counts = ssmcounts)
plotSsm2ss <- function(fit, ssm_counts, linear = TRUE,
                       plot_options = list(colors = NULL, pers_theme = NULL)){

    if(is.null(plot_options$colors)){
        plot_options$colors <- c("#71BD5E", "#B590C7")
    }

    # Create a data frame with the posterior state
    post_stt <- data.frame(mu_est = c(fit$mean$mu_t, fit$mean$mu_wt),
                           mu_lb = c(fit$q2.5$mu_t, fit$q2.5$mu_wt),
                           mu_ub = c(fit$q97.5$mu_t, fit$q97.5$mu_wt),
                           year = rep(unique(ssm_counts$year), 2),
                           season = rep(unique(ssm_counts$season_id), each = nrow(ssm_counts)/2),
                           count = c(log(ssm_counts[ssm_counts$season_id==1, "count", drop = T] + 0.1),
                                     log(ssm_counts[ssm_counts$season_id==2, "count", drop = T] + 0.1)))

    if(linear){
        post_stt <- post_stt %>%
            dplyr::mutate(dplyr::across(.cols = -c(year, season),
                                        .fns = ~exp(.x) - 0.1))
        abund_label <- "Abundance"

        # Cut axis when values are larger than 10 times the max count
        ylims <- c(0, max(post_stt$count, na.rm = T) * 10)

    } else {

        abund_label <- "log abundance"
        ylims <- c(NA, NA)

    }


    # Abundance by season -----------------------------------------------------

    stt_plot <- post_stt %>%
        tidyr::pivot_longer(cols = c(mu_est, mu_lb, mu_ub),
                            names_to = "quantile") %>%
        ggplot2::ggplot() +
        ggplot2::geom_path(aes(x = year, y = value, linetype = quantile)) +
        ggplot2::geom_point(aes(x = year, y = count, col = factor(season)), show.legend = FALSE) +
        ggplot2::scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
        ggplot2::scale_colour_manual(values = plot_options$colors) +
        ggplot2::coord_cartesian(ylim = ylims) +
        ggplot2::facet_wrap("season", nrow = 2,
                            labeller = ggplot2::labeller(season = c("1" = "Summer", "2" = "Winter"))) +
        ggplot2::xlab("Year") + ggplot2::ylab(abund_label) +
        plot_options$pers_theme


    # Trend plot --------------------------------------------------------------

    # Create a data frame with the posterior trend
    post_trd <- data.frame(beta_est = fit$mean$beta,
                           beta_lb = fit$q2.5$beta,
                           beta_ub = fit$q97.5$beta,
                           prop_est = fit$mean$lambda,
                           prop_lb = fit$q2.5$lambda,
                           prop_ub = fit$q97.5$lambda,
                           year = unique(ssm_counts$year))

    if(linear){
        post_trd <- post_trd %>%
            dplyr::mutate(dplyr::across(.cols = -year,
                                        .fns = ~exp(.x)))
        rate_label <- "Growth rate"
        ratio_label <- "W/S ratio"

    } else {

        rate_label <- "log growth rate"
        ratio_label <- "log W/S ratio"

    }

    # Plot trend
    trd_plot <- post_trd %>%
        dplyr::select(-dplyr::starts_with("prop")) %>%
        tidyr::pivot_longer(cols = -year,
                            names_to = "quantile") %>%
        ggplot2::ggplot() +
        ggplot2::geom_path(aes(x = year, y = value, linetype = quantile)) +
        ggplot2::scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
        ggplot2::scale_colour_manual(values = plot_options$colors) +
        ggplot2::xlab("Year") + ggplot2::ylab(rate_label) +
        plot_options$pers_theme


    # Winter/summer ratio -----------------------------------------------------

    # Plot proportion summer/winter
    prop_plot <- post_trd %>%
        dplyr::select(-dplyr::starts_with("beta")) %>%
        tidyr::pivot_longer(cols = -year,
                            names_to = "quantile") %>%
        ggplot2::ggplot() +
        ggplot2::geom_path(aes(x = year, y = value, linetype = quantile)) +
        ggplot2::scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
        ggplot2::scale_colour_manual(values = plot_options$colors) +
        ggplot2::xlab("Year") + ggplot2::ylab(ratio_label) +
        plot_options$pers_theme


    # Save --------------------------------------------------------------------

    # Title
    plottitle <- ifelse(unique(ssm_counts$spp) == "multi",
                        "Multiple species",
                        unique(ssm_counts$spp))

    # To prevent opening multiple devices
    pfile <- tempfile()
    grDevices::png(pfile)
    p <- gridExtra::grid.arrange(stt_plot, trd_plot, prop_plot,
                                 nrow = 3, heights = c(2/4, 1/4, 1/4),
                                 top = plottitle)
    grDevices::dev.off()
    unlink(pfile)

    return(list(plot = p,
                data = list(post_stt, post_trd)))


}
