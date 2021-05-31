#' Plot state series of a CWAC state-space model
#'
#' @param fit A JAGS state-space model fitted to CWAC data
#' @param ssm_counts A data frame with the count data use to fit the state-space model
#' @param dyn Whether the fitted long-term trend is fixed or dynamic (NOT USED AT PRESENT).
#'
#' @return A plot with summer and winter fitted states, as well as the long-term trend
#' @importFrom ggplot2 aes
#' @export
#'
#' @examples
#' counts <- barberspan
#' ssmcounts <- prepSsmData(counts, species = NULL)
#' fit_fxd <- fitCwacSsm2ss(ssmcounts,
#'                          mod_file = "analysis/models/cwac_ssm_2ss_fxd.jags",
#'                          param = c("beta", "sig.w", "sig.eps", "sig.alpha",
#'                                    "sig.e", "mu_t", "mu_wt"))
#' plotSsm(fit = fit_fxd, ssm_counts = ssmcounts)
#'
plotSsm <- function(fit, ssm_counts, dyn = FALSE){

    # Create a data frame with the posterior state
    post_stt <- data.frame(mu_est = c(fit$mean$mu_t, fit$mean$mu_wt),
                                 mu_lb = c(fit$q2.5$mu_t, fit$q2.5$mu_wt),
                                 mu_ub = c(fit$q97.5$mu_t, fit$q97.5$mu_wt),
                                 year = rep(unique(ssm_counts$year), 2),
                                 season = rep(unique(ssm_counts$season_id), each = nrow(ssm_counts)/2),
                                 count = c(log(ssm_counts[ssm_counts$season_id==1, "count", drop = T]),
                                           log(ssm_counts[ssm_counts$season_id==2, "count", drop = T]))) %>%
        dplyr::group_by(year, season) %>%
        dplyr::mutate(seas_id = dplyr::cur_group_id()) %>%
        dplyr::ungroup()

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
                                 trd_ub = fit$q97.5$beta)

        trd_plot <- data.frame(beta = fit$sims.list$beta) %>%
            ggplot2::ggplot() +
            ggplot2::geom_density(aes(x = beta), fill = "yellow") +
            ggplot2::geom_pointrange(data = post_trd, aes(x = trd_est, xmin = trd_lb, xmax = trd_ub, y = 0)) +
            ggplot2::xlab("Trend (log growth-rate)")

    } else {

        # Create a data frame with the posterior trend
        post_trd <- data.frame(beta_est = fit$mean$beta,
                               beta_lb = fit$q2.5$beta,
                               beta_ub = fit$q97.5$beta,
                               year = unique(ssm_counts$year))

        # Plot separated by season
        trd_plot <- post_trd %>%
            tidyr::pivot_longer(cols = -year,
                                names_to = "quantile") %>%
            ggplot2::ggplot() +
            ggplot2::geom_path(aes(x = year, y = value, linetype = quantile)) +
            ggplot2::scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            ggplot2::xlab("Year") + ggplot2::ylab("log growth-rate")
    }

    # Title
    plottitle <- ifelse(unique(ssm_counts$spp) == "multi",
                        "Multiple species",
                        unique(ssm_counts$spp))

    gridExtra::grid.arrange(stt_plot, trd_plot,
                            nrow = 2, heights = c(2/3, 1/3),
                            top = plottitle)

    }
