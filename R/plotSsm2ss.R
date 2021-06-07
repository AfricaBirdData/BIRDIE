#' Plot state series of a CWAC state-space model with two seasons
#'
#' @param fit A JAGS state-space model fitted to CWAC data
#' @param ssm_counts A data frame with the count data use to fit the state-space model
#' @param dyn Whether the fitted long-term trend is fixed or dynamic (NOT USED AT PRESENT).
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
#' plotSsm2ss(fit = fit, ssm_counts = ssmcounts, dyn = TRUE)
plotSsm2ss <- function(fit, ssm_counts, dyn = FALSE){

    # Create a data frame with the posterior state
    post_stt <- data.frame(mu_est = c(fit$mean$mu_t, fit$mean$mu_wt),
                           mu_lb = c(fit$q2.5$mu_t, fit$q2.5$mu_wt),
                           mu_ub = c(fit$q97.5$mu_t, fit$q97.5$mu_wt),
                           year = rep(unique(ssm_counts$year), 2),
                           season = rep(unique(ssm_counts$season_id), each = nrow(ssm_counts)/2),
                           count = c(log(ssm_counts[ssm_counts$season_id==1, "count", drop = T]),
                                     log(ssm_counts[ssm_counts$season_id==2, "count", drop = T])))

    # Plot separated by season
    stt_plot <- post_stt %>%
        tidyr::pivot_longer(cols = c(mu_est, mu_lb, mu_ub),
                     names_to = "quantile") %>%
        ggplot2::ggplot() +
        ggplot2::geom_path(aes(x = year, y = value, linetype = quantile)) +
        ggplot2::geom_point(aes(x = year, y = count, col = factor(season)), show.legend = FALSE) +
        # ggplot2::scale_colour_discrete(name = "Season", labels = c("Summer", "Winter")) +
        ggplot2::scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
                                       #labels = c("Mean", "2.5%", "97.5%")) +
        ggplot2::facet_wrap("season", nrow = 2,
                            labeller = ggplot2::labeller(season = c("1" = "Summer", "2" = "Winter"))) +
        ggplot2::xlab("Year") + ggplot2::ylab("log-abundance")

    # Create a plot for the trend
    if(dyn == FALSE){

        # Create a data frame with the posterior trend
        post_trd <- data.frame(trd_est = fit$mean$beta,
                               trd_lb = fit$q2.5$beta,
                               trd_ub = fit$q97.5$beta,
                               prop_est = fit$mean$lambda,
                               prop_lb = fit$q2.5$lambda,
                               prop_ub = fit$q97.5$lambda)

        trd_plot <- data.frame(beta = fit$sims.list$beta) %>%
            ggplot2::ggplot() +
            ggplot2::geom_density(aes(x = beta), fill = "yellow") +
            ggplot2::geom_pointrange(data = post_trd, aes(x = trd_est, xmin = trd_lb, xmax = trd_ub, y = 0)) +
            ggplot2::xlab("Trend (log growth-rate)")

        prop_plot <- data.frame(lambda = fit$sims.list$lambda) %>%
            ggplot2::ggplot() +
            ggplot2::geom_density(aes(x = lambda), fill = "yellow") +
            ggplot2::geom_pointrange(data = post_trd, aes(x = trd_est, xmin = trd_lb, xmax = trd_ub, y = 0)) +
            ggplot2::xlab("Log propotion summer/winter")

    } else {

        # Create a data frame with the posterior trend
        post_trd <- data.frame(beta_est = fit$mean$beta,
                               beta_lb = fit$q2.5$beta,
                               beta_ub = fit$q97.5$beta,
                               prop_est = fit$mean$lambda,
                               prop_lb = fit$q2.5$lambda,
                               prop_ub = fit$q97.5$lambda,
                               year = unique(ssm_counts$year))

        # Plot trend
        trd_plot <- post_trd %>%
            dplyr::select(-dplyr::starts_with("prop")) %>%
            tidyr::pivot_longer(cols = -year,
                                names_to = "quantile") %>%
            ggplot2::ggplot() +
            ggplot2::geom_path(aes(x = year, y = value, linetype = quantile)) +
            ggplot2::scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            ggplot2::xlab("Year") + ggplot2::ylab("log growth-rate")

        # Plot proportion summer/winter
        prop_plot <- post_trd %>%
            dplyr::select(-dplyr::starts_with("beta")) %>%
            tidyr::pivot_longer(cols = -year,
                                names_to = "quantile") %>%
            ggplot2::ggplot() +
            ggplot2::geom_path(aes(x = year, y = value, linetype = quantile)) +
            ggplot2::scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            ggplot2::xlab("Year") + ggplot2::ylab("Log proportion")
    }

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
