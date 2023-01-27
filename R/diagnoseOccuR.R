#' Run diagnostics for occuR model
#'
#' @description This function runs basic diagnostic checks for occuR models.
#' Routines are based on those in spOccupancy and perform posterior predictive
#' goodness of fit checks see \code{\link[spOccupancy]{ppcOcc}}. For occuR predictions
#' are obtained via bootstrap.
#' @inheritParams ppl_run_pipe_dst1
#' @param fit An occuR model fit see \code{\link[occuR]{fit_occu}}
#' @param year_sel The year the model was run for.
#'
#' @return The function will create a list with several objects. The most important are:
#' - fit.y: For each MCMC sample, chi-squared statistic for total number of detections
#' for all sites in data wrt expected model expectation
#' - fit.y.rep: For each MCMC sample, chi-squared statistic for total number of detections
#' for all sites in simulated data wrt expected model expectation
#' - fit.y.group.quants: For each site, chi-squared statistic for total number of detections
#' in data wrt expected model expectation (posterior distribution quantiles)
#' - fit.y.rep.group.quants: For each site, chi-squared statistic for total number of detections
#' in simulated data wrt expected model expectation (posterior distribution quantiles)
#' - y.summ.per.site: Data frame with number of detections per site obtained from
#' data simulated from the posterior distribution. Also the observed number of detections
#' in the data.
#'
#' @export
#'
#' @examples
diagnoseOccuR <- function(fit, sp_code, config, year_sel){

    # Check convergence
    fails <- fit$fit$convergence

    if(fails != 0){
        conv_file <- file.path(config$out_dir, paste0("reports/no_converge_occu_fit_", config$package, "_", year_sel, "_", sp_code, ".txt"))
        sink(conv_file)
        print(paste("non-convergence", fit$fit$message))
        sink()
        message(paste("non-convergence", fit$fit$message)) # to console
    }

    # Simulate detections/non-detection data from model
    post_sims <- simDetOccuR(fit)

    # Check goodness-of-fit (Bayesian p-value from posterior predictive check)
    ppc_out <- gofOccuR(fit, post_sims, fit_stat = 'chi-squared', group = 1)

    bayes_p <- mean(ppc_out$fit.y.rep > ppc_out$fit.y)

    if(bayes_p < 0.05){
        sink(file.path(config$out_dir, paste0("reports/gof_occu_fit_", config$package, "_", year_sel, "_", sp_code, ".txt")))
        print(paste("Bayesian p-value =", round(bayes_p, 2)))
        sink()
        message(paste("Bayesian p-value =", round(bayes_p, 2)))
    }

    # Posterior coverage
    y_sims_per_site <- as.data.frame(post_sims$y.rep.samples) %>%
        dplyr::mutate(site_id = attr(post_sims$y.rep.samples, "site_id")) %>%
        dplyr::group_by(site_id) %>%
        dplyr::summarise(dplyr::across(.cols = dplyr::starts_with("V"), .fns = ~sum(.x))) %>%
        dplyr::ungroup() %>%
        dplyr::select(-1) %>%
        as.matrix()

    y_sims_per_site_quants <- apply(y_sims_per_site, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))

    y_grouped <- data.frame(y = fit$mats$y) %>%
        dplyr::mutate(site_id = attr(fit$mats$X_p, "site")) %>%
        dplyr::group_by(site_id) %>%
        dplyr::summarise(n = sum(y)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-1) %>%
        as.matrix()

    pst_df <- as.data.frame(t(y_sims_per_site_quants)) %>%
        stats::setNames(paste0("q", c(0.025, 0.25, 0.5, 0.75, 0.975)*100)) %>%
        dplyr::mutate(obs = y_grouped,
                      site_id = names(y_grouped))

    # Save per site info
    ppc_out$y.summ.per.site <- pst_df

    return(ppc_out)

}

