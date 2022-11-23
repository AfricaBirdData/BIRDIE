#' Prepare spOccupancy data for model fitting
#'
#' @inheritParams ppl_run_pipe_dst1
#' @param spatial Logical, indicating whether spatial random effects should be
#' included in the model. Defaults to FALSE.
#'
#' @return
#' @export
#'
#' @examples
prepSpOccuData <- function(sp_code, year, config, spatial = FALSE, ...){

    varargs <- list(...)


    # Load data ---------------------------------------------------------------

    # File names
    visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
    sitefile <- file.path(config$out_dir, paste0("occu_site_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

    # Read in site and visit data
    site_data <- utils::read.csv(sitefile)
    visit_data <- utils::read.csv(visitfile)
    det_data <- utils::read.csv(detfile)

    # Stop if there are no detections
    if(!1 %in% unique(det_data$obs)){
        warning(paste("No detection of species", sp_code))
        return(1)
    }

    # Or species detected in too few Pentads
    n_pentads <- det_data %>%
        dplyr::count(Pentad, obs) %>%
        dplyr::filter(obs == 1) %>%
        nrow()

    if(n_pentads < 5){
        warning(paste("Species", sp_code, "detected in less than 5 pentads"))
        return(2)
    }


    # Prepare spOccupancy data list -------------------------------------------

    # Add detection info to visit data
    visit_data <- visit_data %>%
        dplyr::left_join(det_data,
                         by = c("CardNo", "StartDate", "year", "Pentad")) %>%
        dplyr::mutate(Spp = ifelse(obs == 0, "-", 1))

    # Return list for spatial occupancy model
    if(spatial){
        # Add coordinates to the data
        spocc_data <- ABAP::abapToSpOcc_multi(visit_data,
                                          pentads = sa_pentads %>%
                                              dplyr::filter(Name %in% unique(site_data$Name)),
                                          seasons = "year")
    } else {
        spocc_data <- ABAP::abapToSpOcc_multi(visit_data, seasons = "year")
    }


    # Add covariates to spOccupancy object -----------------------------------

    # Keep only those sites that appear in visits
    site_data <- site_data %>%
        dplyr::filter(Pentad %in% unique(visit_data$Pentad))

    # Select covariates and add to data list
    spocc_data <- spocc_data %>%
        ABAP::addEEtoSpOcc_multi(
            ee_data = site_data %>%
                dplyr::select(pentad = Pentad, year, dist_coast, prcp, tdiff, ndvi,
                              watext, watrec, elev, log_dist_coast, log_watext),
            type = "occ", seasons = "year")

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "occ", scale_vars = names(spocc_data$occ.covs)[-1])

    # Add detection covariates
    spocc_data <- spocc_data %>%
        addSpOccDetCovt(visit_data %>%
                            dplyr::mutate(StartMonth = lubridate::month(StartDate)) %>%
                            dplyr::select(Pentad, StartDate, StartMonth, year),
                        seasons = "year")


    # Add observer ID as covariate
    spocc_data <- spocc_data %>%
        addSpOccDetCovt(visit_data %>%
                            dplyr::mutate(obs_id = as.numeric(ObserverNo)) %>%
                            dplyr::select(Pentad, StartDate, obs_id, year),
                        seasons = "year")

    # Add temperature, precipitation, CWAC site presence and pentad as detection covariates
    spocc_data <- spocc_data %>%
        ABAP::addEEtoSpOcc_multi(
            ee_data = site_data %>%
                dplyr::group_by(Pentad) %>%
                dplyr::mutate(site_id = dplyr::cur_group_id()) %>%
                dplyr::ungroup() %>%
                dplyr::select(pentad = Pentad, year, prcp, tdiff, cwac, site_id),
            type = "det", seasons = "year")

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "det", scale_vars = c("prcp", "tdiff"))

    # # Create a cyclic month variable
    # spocc_data$det.covs$month_sin <- sin(2*pi*spocc_data$det.covs$month/12)
    # spocc_data$det.covs$month_cos <- cos(2*pi*spocc_data$det.covs$month/12)


    return(spocc_data)

}
