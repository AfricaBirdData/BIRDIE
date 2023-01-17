#' Prepare occuR data for model fitting
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
prepOccuRData <- function(site_data, visit_data, config, spatial = FALSE, sp_sites = NULL){


    # Prepare occuR data list -------------------------------------------------

    # Return list for spatial occupancy model
    if(spatial){
        # Add coordinates to the data
        occur_data <- ABAP::abapToOccuR(visit_data,
                                        occasion = "year",
                                        pentads = sp_sites %>%
                                            dplyr::filter(Name %in% unique(site_data$Name)))
    } else {
        occur_data <- ABAP::abapToOccuR(visit_data,
                                        occasion = "year")
    }


    # Add covariates to occuR object ------------------------------------------

    # Select occupancy covariates and add to data list
    tt_occ <- stats::terms(stats::reformulate(config$occ_mod))

    occ_vars <- attr(tt_occ, "term.labels")
    occ_vars <- gsub(".* \\| ", "", occ_vars)

    occ_cov_sel <- site_data %>%
        dplyr::select(pentad = Pentad, year, dplyr::all_of(occ_vars))

    occur_data$site <- occur_data$site %>%
        dplyr::left_join(occ_cov_sel, by = c("pentad", "year"))

    # Scale covariates
    occur_data$site <- occur_data$site %>%
        dplyr::mutate(dplyr::across(-c(pentad, year, site, occasion),
                                    ~ scale(.x)))

    # Add detection covariates
    tt_det <- stats::terms(stats::reformulate(config$det_mod))

    det_vars <- attr(tt_det, "term.labels")
    det_vars <- gsub(".* \\| ", "", det_vars)

    det_cov_sel <- visit_data %>%
        dplyr::arrange(Pentad, StartDate) %>%  # This is how abapToOccuR orders rows
        dplyr::select(StartDate, dplyr::all_of(det_vars))

    if(nrow(det_cov_sel) == nrow(occur_data$visit)){
        occur_data$visit <- occur_data$visit %>%
            dplyr::select(pentad, year, site, occasion, visit, obs) %>%
            dplyr::bind_cols(det_cov_sel)
    } else {
        stop("Occupancy visit data and covariate data have different number of rows")
    }

    # Scale covariates
    occur_data$visit <- occur_data$visit %>%
        dplyr::mutate(dplyr::across(c(log_hours, prcp, tdiff),
                                    ~ scale(.x)))

    # # Create a cyclic month variable
    # spocc_data$det.covs$month_sin <- sin(2*pi*spocc_data$det.covs$month/12)
    # spocc_data$det.covs$month_cos <- cos(2*pi*spocc_data$det.covs$month/12)


    return(occur_data)

}
