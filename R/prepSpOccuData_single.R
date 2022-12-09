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
#'
#' @return
#' @export
#'
#' @examples
prepSpOccuData_single <- function(site_data, visit_data, config, spatial = FALSE){


    # Prepare spOccupancy data list -------------------------------------------

    # Return list for spatial occupancy model
    if(spatial){
        # Add coordinates to the data
        spocc_data <- ABAP::abapToSpOcc_single(visit_data,
                                          pentads = sa_pentads %>%
                                              dplyr::filter(Name %in% unique(site_data$Name)))
    } else {
        spocc_data <- ABAP::abapToSpOcc_single(visit_data)
    }


    # Add covariates to spOccupancy object -----------------------------------

    # Keep only those sites that appear in visits
    site_data <- site_data %>%
        dplyr::filter(Pentad %in% unique(visit_data$Pentad))

    # Select covariates and add to data list
    occ_cov_sel <- site_data %>%
        dplyr::select(pentad = Pentad, dist_coast, prcp, tdiff, ndvi,
                      watext, watrec, elev, log_dist_coast, log_watext)

    spocc_data <- spocc_data %>%
        ABAP::addEEtoSpOcc_single(ee_data = occ_cov_sel)

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "occ", scale_vars = names(occ_cov_sel)[-1])

    # Add detection covariates
    spocc_data <- spocc_data %>%
        addSpOccDetCovt(visit_data %>%
                            dplyr::mutate(StartMonth = lubridate::month(StartDate)) %>%
                            dplyr::select(Pentad, StartDate, StartMonth))


    # Add observer ID as covariate
    spocc_data <- spocc_data %>%
        addSpOccDetCovt(visit_data %>%
                            dplyr::mutate(obs_id = as.numeric(ObserverNo)) %>%
                            dplyr::select(Pentad, StartDate, obs_id))

    # Add temperature, precipitation, CWAC site presence and pentad as detection covariates
    spocc_data <- spocc_data %>%
        addSpOccDetCovt(visit_data %>%
                            dplyr::group_by(Pentad) %>%
                            dplyr::mutate(site_id = dplyr::cur_group_id(),
                                          tdiff = tmmx - tmmn) %>%
                            dplyr::ungroup() %>%
                            dplyr::select(Pentad, StartDate, prcp, tdiff, cwac, site_id))

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "det", scale_vars = c("prcp", "tdiff"))

    # # Create a cyclic month variable
    # spocc_data$det.covs$month_sin <- sin(2*pi*spocc_data$det.covs$month/12)
    # spocc_data$det.covs$month_cos <- cos(2*pi*spocc_data$det.covs$month/12)


    return(spocc_data)

}
