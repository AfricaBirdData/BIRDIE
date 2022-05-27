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
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' counts <- barberspan
#' ssmcounts <- prepSsmData(counts, species = NULL)
#' fit <- fitCwacSsm(ssmcounts, mod_file = "mymodel.jags",
#' param = c("beta", "lambda", "sig.zeta",
#' "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))
#' plotSsm2ss(fit = fit, ssm_counts = ssmcounts)
#' }
plotSsm2ss <- function(fit, ssm_counts, linear = TRUE,
                       plot_options = list(colors = NULL, pers_theme = NULL)){

    if(is.null(plot_options$colors)){
        plot_options$colors <- c("#71BD5E", "#B590C7")
    }

    # Create a data frame with the posterior state
    post_stt <- data.frame(mu_est = fit$mean$mu_t,
                           mu_lb = fit$q2.5$mu_t,
                           mu_ub = fit$q97.5$mu_t,
                           year = ssm_counts$year,
                           season = ssm_counts$Season,
                           site = ssm_counts$LocationCode,
                           count = log(ssm_counts$count + 0.1)) %>%
        dplyr::filter(season != "O") %>%
        dplyr::mutate(season = ifelse(season == "S", 1, 2))

    # There might be counts in more than one day per season and year, so
    # we plot the mean for the season
    post_stt <- post_stt %>%
        dplyr::group_by(year, season, site) %>%
        dplyr::summarize(count = mean(count),
                         mu_est = mean(mu_est),
                         mu_lb = mean(mu_lb),
                         mu_ub = mean(mu_ub)) %>%
        dplyr::ungroup()

    if(linear){
        post_stt <- post_stt %>%
            dplyr::mutate(dplyr::across(.cols = -c(year, season, site),
                                        .fns = ~exp(.x) - 0.1))
        abund_label <- "Abundance"

        # Cut axis when values are larger than 10 times the max count
        # ylims <- c(0, max(post_stt$count, na.rm = T) * 10)

    } else {

        abund_label <- "log abundance"
        ylims <- c(NA, NA)

    }


    site_plot_data <- vector("list", length(unique(post_stt$site)))

    for(i in seq_along(site_plot_data)){

        # Abundance by season -----------------------------------------------------

        stt_plot <- post_stt %>%
            dplyr::filter(site == unique(post_stt$site)[i]) %>%
            tidyr::pivot_longer(cols = c(mu_est, mu_lb, mu_ub),
                                names_to = "quantile") %>%
            ggplot() +
            geom_path(aes(x = year, y = value, linetype = quantile)) +
            geom_point(aes(x = year, y = count, col = factor(season)), show.legend = FALSE) +
            scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            scale_colour_manual(values = plot_options$colors) +
            # coord_cartesian(ylim = ylims) +
            facet_wrap("season", ncol = 2,
                       labeller = labeller(season = c("1" = "Summer", "2" = "Winter"))) +
            xlab("Year") + ylab(abund_label) +
            plot_options$pers_theme


        # Trend plot --------------------------------------------------------------

        # Create a data frame with the posterior trend
        post_trd <- data.frame(beta_est = fit$mean$beta[i,],
                               beta_lb = fit$q2.5$beta[i,],
                               beta_ub = fit$q97.5$beta[i,],
                               prop_est = fit$mean$lambda[i,],
                               prop_lb = fit$q2.5$lambda[i,],
                               prop_ub = fit$q97.5$lambda[i,],
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
            ggplot() +
            geom_path(aes(x = year, y = value, linetype = quantile)) +
            scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            scale_colour_manual(values = plot_options$colors) +
            xlab("Year") + ylab(rate_label) +
            plot_options$pers_theme


        # Winter/summer ratio -----------------------------------------------------

        # Plot proportion summer/winter
        prop_plot <- post_trd %>%
            dplyr::select(-dplyr::starts_with("beta")) %>%
            tidyr::pivot_longer(cols = -year,
                                names_to = "quantile") %>%
            ggplot() +
            geom_path(aes(x = year, y = value, linetype = quantile)) +
            scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            scale_colour_manual(values = plot_options$colors) +
            xlab("Year") + ylab(ratio_label) +
            plot_options$pers_theme


        # Save --------------------------------------------------------------------

        # Title
        plottitle <- paste(ifelse(unique(ssm_counts$spp) == "multi",
                                  "Multiple species",
                                  unique(ssm_counts$spp)),
                           "at site",
                           unique(ssm_counts$LocationCode)[i])

        # To prevent opening multiple devices
        pfile <- tempfile()
        grDevices::png(pfile)
        p <- gridExtra::grid.arrange(stt_plot, trd_plot, prop_plot,
                                     layout_matrix = matrix(c(1,2,1,3), nrow = 2),
                                     top = plottitle)
        grDevices::dev.off()
        unlink(pfile)

        site_plot_data[[i]] <- list(plot = p,
                                    data = list(post_stt %>%
                                                    dplyr::filter(site == unique(post_stt$site)[i]), post_trd))

    }

    names(site_plot_data) <- unique(ssm_counts$LocationCode)

    return(site_plot_data)


}
