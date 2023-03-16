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

    # Data summary
    dat_summ <- ssm_counts %>%
        dplyr::mutate(season_est = as.integer(fit$mean$summer < 0.5) + 1L) %>%
        dplyr::group_by(site_id, year, season_est) %>%
        dplyr::summarise(count = log(mean(count+1, na.rm = TRUE))) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = season_est, values_from = count) %>%
        dplyr::rename(count_s = "1",
                      count_w = "2")


    # Create a data frame with the posterior state
    post_stt <- data.frame(year = dat_summ$year,
                           site_id = dat_summ$site_id,
                           stt_s_est = as.vector(t(fit$mean$stt_s)),
                           stt_s_lb = as.vector(t(fit$q2.5$stt_s)),
                           stt_s_ub = as.vector(t(fit$q97.5$stt_s)),
                           stt_w_est = as.vector(t(fit$mean$stt_w)),
                           stt_w_lb = as.vector(t(fit$q2.5$stt_w)),
                           stt_w_ub = as.vector(t(fit$q97.5$stt_w)))

    # Pivot seasons
    dat_summ <- dat_summ %>%
        tidyr::pivot_longer(cols = c(count_s, count_w),
                            names_to = "season", values_to = "count") %>%
        dplyr::mutate(season = ifelse(grepl("_s", season), "summer", "winter"))

    post_stt <- post_stt %>%
        tidyr::pivot_longer(cols = -c(year, site_id),
                            names_to = "quantile", values_to = "value") %>%
        dplyr::mutate(season = ifelse(grepl("_s", quantile), "summer", "winter")) %>%
        dplyr::mutate(quantile = dplyr::case_when(grepl("_est", quantile) ~ "est",
                                                  grepl("_ub", quantile) ~ "ub",
                                                  grepl("_lb", quantile) ~ "lb"))

    # Add counts and location code
    post_stt <- post_stt %>%
        dplyr::left_join(dplyr::distinct(ssm_counts[,c("site_id", "LocationCode")]),
                         by = "site_id") %>%
        dplyr::rename(loc_code = LocationCode) %>%
        dplyr::left_join(dat_summ, by = c("site_id", "year", "season"))



    if(linear){
        post_stt <- post_stt %>%
            dplyr::mutate(value = exp(value) - 1,
                          count = exp(count) - 1) %>%
            dplyr::mutate(value = ifelse(value < 0, 0, value))

        abund_label <- "Abundance"

        # Cut axis when values are larger than 10 times the max count
        # ylims <- c(0, max(post_stt$count, na.rm = T) * 10)

    } else {

        abund_label <- "log abundance"
        ylims <- c(NA, NA)

    }


    site_plot_data <- vector("list", dplyr::n_distinct(post_stt$loc_code))


    for(i in seq_along(site_plot_data)){

        loc_code_sel <- post_stt %>%
            dplyr::filter(site_id == i) %>%
            dplyr::distinct(loc_code) %>%
            dplyr::pull(loc_code)

        # Abundance by season -----------------------------------------------------

        stt_plot <- post_stt %>%
            dplyr::filter(site_id == i) %>%
            ggplot() +
            geom_path(aes(x = year, y = value, linetype = quantile)) +
            geom_point(aes(x = year, y = count, col = season), show.legend = FALSE) +
            scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            scale_colour_manual(values = plot_options$colors) +
            # coord_cartesian(ylim = ylims) +
            facet_wrap("season", ncol = 2) +
            xlab("Year") + ylab(abund_label) +
            plot_options$pers_theme


        # Trend plot --------------------------------------------------------------

        # Create a data frame with the posterior trend
        post_trd <- data.frame(beta_est = c(fit$mean$beta[i,], NA),
                               beta_lb = c(fit$q2.5$beta[i,], NA),
                               beta_ub = c(fit$q97.5$beta[i,], NA),
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
                           loc_code_sel)

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
                                                    dplyr::filter(site_id == i),
                                                post_trd))

        names(site_plot_data)[[i]] <- loc_code_sel

    }

    return(site_plot_data)


}
