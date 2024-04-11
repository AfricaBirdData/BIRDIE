#' Run diagnostics for occupancy model fit
#'
#' @description This function performs basic checks for an occupancy model fit,
#' such as convergence of parameters and posterior predictive checks
#' @inheritParams ppl_run_pipe_dst1
#' @param fit A `spOccupancy` model fit to ABAP data.
#' @param ppc A list with the outputs from \link{diagnoseGofSpOccu}.
#' @param data A dataset with detection/non-detection data. Not used at the moment because
#' spOccupancy (currently used) fit objects contain the data used to fit the model.
#' @param year The year the pipeline is run for
#'
#' @return A data frame with Rhat values for the different parameters estimated
#' by the model is returned. If `ppc` is provided it will save some time because
#' posterior simulations can be time consuming. Sometimes the posterior GOF test
#' and posterior simulations might be available from previous pipeline runs.
#'
#' @export
#'
#' @examples
ppl_diagnose_occu <- function(fit, ppc = NULL, data = NULL, sp_code, year){

    diag_out <- list()

    # Diagnose convergence
    diag_out$rhat <- diagnoseRhatSpOccu(fit, sp_code, year)

    # Add original number of visits to get an idea of data lost to NA covariates
    year_sel <- year
    geevisitfile <- file.path(config$out_dir, paste0("visit_dat_sa_gee_", config$years_ch, ".csv"))
    geevisitdata <- utils::read.csv(geevisitfile)
    diag_out$rhat$nvisits0 <- geevisitdata |>
        dplyr::filter(year == year_sel) |>
        nrow()

    # Diagnose model fit
    if(is.null(ppc)){
        post_sims <- simDetSpOccu(fit)
        diag_out$gof <- diagnoseGofSpOccu(fit, post_sims, fit_stat = "chi-squared", group = 1)
    } else {
        diag_out$gof <- ppc
    }


    return(diag_out)

}
