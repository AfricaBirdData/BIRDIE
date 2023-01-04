#' Prepare spOccupancy data for single-season model fitting
#'
#' @inheritParams ppl_run_pipe_dst1
#' @param site_data A data frame containing information about covariates associated
#' with ABAP pentads for a given year.
#' @param visit_data A data frame with information associated with sampling visits
#' to ABAP pentads in a given year. Detection/non-detection data for the species
#' of interest must also be included in this data frame.
#' @param spatial Logical, indicating whether spatial random effects should be
#' included in the model. Defaults to FALSE.
#' @sp_sites Spatial object with the sites to be used for fitting spatial models
#'
#' @return
#' @export
#'
#' @examples
prepSpOccuData_single <- function(site_data, visit_data, config, spatial = FALSE, sp_sites = NULL){


    # Prepare spOccupancy data list -------------------------------------------

    # Return list for spatial occupancy model
    if(spatial){
        # Add coordinates to the data
        spocc_data <- ABAP::abapToSpOcc_single(visit_data,
                                          pentads = sp_sites %>%
                                              dplyr::filter(Name %in% unique(site_data$Name)))
    } else {
        spocc_data <- ABAP::abapToSpOcc_single(visit_data)
    }


    # Add covariates to spOccupancy object -----------------------------------

    # Select occupancy covariates and add to data list
    tt_occ <- stats::terms(stats::reformulate(config$occ_mod))

    occ_vars <- attr(tt_occ, "term.labels")
    occ_vars <- gsub(".* \\| ", "", occ_vars)

    occ_cov_sel <- site_data %>%
        dplyr::select(pentad = Pentad, dplyr::all_of(occ_vars))

    spocc_data <- spocc_data %>%
        ABAP::addEEtoSpOcc_single(ee_data = occ_cov_sel)

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "occ", scale_vars = names(occ_cov_sel)[-1])

    # Add detection covariates
    tt_det <- stats::terms(stats::reformulate(config$det_mod))

    det_vars <- attr(tt_det, "term.labels")
    det_vars <- gsub(".* \\| ", "", det_vars)

    det_cov_sel <- visit_data %>%
        dplyr::select(Pentad, StartDate, dplyr::all_of(det_vars))

    spocc_data <- addSpOccDetCovt(spocc_data, det_cov_sel)

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "det", scale_vars = c("log_hours", "prcp", "tdiff"))

    # # Create a cyclic month variable
    # spocc_data$det.covs$month_sin <- sin(2*pi*spocc_data$det.covs$month/12)
    # spocc_data$det.covs$month_cos <- cos(2*pi*spocc_data$det.covs$month/12)


    return(spocc_data)

}
