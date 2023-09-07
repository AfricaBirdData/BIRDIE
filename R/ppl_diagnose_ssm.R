#' Run diagnostics for a state-space model fit
#'
#' @description This function performs basic checks for an SSM fit, such as
#' convergence of parameters and posterior predictive checks
#' @inheritParams ppl_run_pipe_abu1
#' @param fit A JAGS state-space model fitted to CWAC data
#' @param counts A dataframe with CWAC counts ready for model fit. See \code{\link{ppl_create_data_ssm}}
#'
#' @return A list with Rhat values and posterior check statistics is returned.
#' At them moment, we obtain Rhat values for all monitored parameters and three
#' posterior check statistics: "Tmean" proportion of posterior simulations with mean
#' greater than that observed in the data (we would like values close to 0.5),
#' "Tsd" proportion of posterior simulations with sd greater than that observed in
#' the data (we would like values close to 0.5), "Tdiff" mean difference between
#' observed data and posterior simulations (we would like values close to 0).
#'
#' @export
#'
#' @examples
ppl_diagnose_ssm <- function(fit, counts, sp_code, config){

    # Diagnose convergence
    diag_df <- diagnoseRhatJagsSsm(fit, sp_code, config)

    # Diagnose model fit
    gof <- diagnoseGofJagsSsm(fit, counts)
    gof <- round(gof, 3)

    diag_df$Tmean <- gof["Tmean"]
    diag_df$Tsd <- gof["Tsd"]
    diag_df$Tdiff <- gof["Tdiff"]

    return(diag_df)

}
