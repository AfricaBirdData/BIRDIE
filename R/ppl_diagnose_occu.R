#' Run diagnostics for occupancy model fit
#'
#' @description This function performs basic checks for an occupancy model fit,
#' such as convergence of parameters and posterior predictive checks
#' @inheritParams ppl_run_pipe_dst1
#' @param fit A `spOccupancy` model fit to ABAP data.
#' @param data A dataset with detection/non-detection data. Not used at the moment because
#' spOccupancy (currently used) fit objects contain the data used to fit the model.
#' @param year The year the pipeline is run for
#'
#' @return A data frame with Rhat values for the different parameters estimated
#' by the model is returned.
#'
#' @export
#'
#' @examples
ppl_diagnose_occu <- function(fit, data = NULL, sp_code, year){

    diag_out <- list()

    # Diagnose convergence
    diag_out$rhat <- diagnoseRhatSpOccu(fit, sp_code, year)

    # Diagnose model fit
    post_sims <- simDetSpOccu(fit)
    diag_out$gof <- diagnoseGofSpOccu(fit, post_sims, fit_stat = "chi-squared", group = 1)

    return(diag_out)

}
