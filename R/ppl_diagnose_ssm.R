#' Run diagnostics for JAGS state-space model fit
#'
#' @description This function performs basic checks for a JAGS SSM fit, such as
#' convergence of parameters and posterior predictive checks
#' @inheritParams ppl_run_pipe_abu1
#'
#' @return A data frame with Rhat values for the different parameters estimated
#' by the model, will be saved to disk.
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
