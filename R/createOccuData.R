#' Create data for fitting an occupancy model
#'
#' @description This function takes generic occupancy site data and, optionally,
#' visit data without species observations, adds detection data for a species in a given
#' year, runs some checks and formats the output.
#' @param sp_code SAFRING_No of the species of interest as extracted from ABAP.
#' Ignored if download set to FALSE.
#' @param years A numeric vector with elements corresponding to the years we
#' want data for. Ignored if download set to FALSE.
#' @param site_data A dataframe with occupancy site data.
#' @param visit_data Optional. A dataframe with occupancy visit data.
#' @param config A list of configuration parameters see \link{configPreambOccu}
#'
#' @return A list containing two data frames: one with site data and one with
#' visit data, if visit data is provided. The second element will be NULL if visit
#' data not provided. Note that if `visit_data` is provided then all the sites without
#' visits will be removed - i.e., both site and visit data frames will have the same
#' set of sites.
#' @export
#'
#' @examples
createOccuData <- function(config, sp_code, years,
                           site_data, visit_data = NULL){

    # Configure site covariates -----------------------------------------------

    # Expand time-varying covariates and remove those rows with NA values
    site_data <- site_data %>%
        dplyr::rename(Pentad = Name) %>%
        dplyr::select(dplyr::all_of(config$fixed_vars),
                      dplyr::ends_with(match = as.character(years))) %>%
        BIRDIE::gatherYearFromVars(vars = setdiff(names(.), config$fixed_vars), sep = "_") %>%
        tidyr::drop_na()   # I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE

    # Create other variables that could be necessary for the model
    site_data <- site_data %>%
        dplyr::arrange(Pentad) %>%
        dplyr::group_by(Pentad) %>%
        dplyr::mutate(tdiff = tmmx - tmmn,
                      log_dist_coast = log(dist_coast + 1),
                      log_watext = log(watext + 1),
                      log_watrec = log(watrec + 1),
                      log_wetext = log(wetext_2018 + 1),
                      wetcon = wetcon_2018,
                      site_id = dplyr::cur_group_id()) %>%
        dplyr::ungroup()

    # Include interaction terms if there are any
    tt_occu <- stats::terms(stats::reformulate(config$occ_mod))

    intrcs <- attr(tt_occu, "term.labels")[attr(tt_occu, "order") > 1]
    intrcs_vars <- strsplit(intrcs, ":")

    for(i in seq_along(intrcs)){
        site_data <- site_data %>%
            dplyr::mutate(!!intrcs[i] := !!rlang::sym(intrcs_vars[[i]][1]) * !!rlang::sym(intrcs_vars[[i]][2]))
    }


    if(!is.null(visit_data)){

        # Configure visit covariates ----------------------------------------------

        # Subset visits that correspond to the years of interest
        visit_data <- visit_data %>%
            dplyr::filter(year %in% years) %>%
            tidyr::drop_na()   # I'M REMOVING VISITS WITH NA DATA! MAKE SURE THIS MAKES SENSE

        # Create additional detection variables
        visit_data <- visit_data  %>%
            dplyr::mutate(StartDate = as.Date(StartDate))%>%
            dplyr::mutate(StartMonth = lubridate::month(StartDate),
                          obs_id = as.numeric(ObserverNo),
                          tdiff = tmmx - tmmn,
                          hours = TotalHours,
                          log_hours = log(hours + 1)) %>%
            dplyr::left_join(site_data %>%
                                 dplyr::select(Pentad, site_id) %>%
                                 dplyr::distinct(),
                             by = "Pentad")

        # Include interaction terms if there are any
        tt_det <- stats::terms(stats::reformulate(config$det_mod))

        intrcs <- attr(tt_det, "term.labels")[attr(tt_det, "order") > 1]
        intrcs_vars <- strsplit(intrcs, ":")

        for(i in seq_along(intrcs)){
            visit_data <- visit_data %>%
                dplyr::mutate(!!intrcs[i] := !!rlang::sym(intrcs_vars[[i]][1]) * !!rlang::sym(intrcs_vars[[i]][2]))
        }

    }

    # Create other covariates -------------------------------------------------

    message("Creating additional covariates...")

    # Load South African pentads
    sa_pentads <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa")

    # Create a detection covariate for CWAC sites
    cwac_sites <- CWAC::listCwacSites(.region_type = "country",
                                      .region = "South Africa") %>%
        sf::st_as_sf(coords = c("X", "Y"), crs = 4326, remove = FALSE) %>%
        suppressWarnings()

    # Make sites a spatial object
    site_data <- site_data %>%
        dplyr::left_join(sa_pentads %>%
                             dplyr::select(Name),
                         by = c("Pentad" = "Name")) %>%
        sf::st_sf()

    cwac_intsc <- sf::st_intersects(site_data, cwac_sites)
    site_data <- site_data %>%
        dplyr::mutate(cwac = as.integer(sapply(cwac_intsc, length) != 0)) %>%
        sf::st_drop_geometry()

    if(!is.null(visit_data)){

        # Add cwac to visit data
        visit_data <- visit_data %>%
            dplyr::left_join(site_data %>%
                                 dplyr::select(Pentad, cwac) %>%
                                 dplyr::distinct(),
                             by = "Pentad")

        det_vars <- attr(tt_det, "term.labels")
        det_vars <- gsub(".* \\| ", "", det_vars)

        # We need to keep some extra variables of the visit data for other functions to work down the line
        visit_data <- visit_data %>%
            dplyr::select(CardNo, Pentad, site_id, StartDate, year, TotalHours, dplyr::all_of(det_vars))

    }

    # Keep only those variables that will be used in the models
    occu_vars <- attr(tt_occu, "term.labels")
    occu_vars <- gsub(".* \\| ", "", occu_vars)

    site_data <- site_data %>%
        dplyr::select(Pentad, site_id, year, dplyr::all_of(occu_vars))

    # Clean up data -----------------------------------------------------------

    if(!is.null(visit_data)){

        # Keep only pentads that appear in visit data
        site_data <- site_data %>%
            dplyr::filter(Pentad %in% unique(visit_data$Pentad))

        # Check that all pentads-years in visit data are also in site data
        miss_pentads <- visit_data %>%
            dplyr::as_tibble() %>%
            dplyr::distinct(Pentad, year) %>%
            dplyr::anti_join(
                site_data %>%
                    dplyr::as_tibble() %>%
                    dplyr::distinct(Pentad, year),
                by = c("Pentad", "year"))

        if(nrow(miss_pentads) > 0){
            log_name <- file.path(config$out_dir, "reports", paste0("pentads_in_visit_not_site_", config$years_ch, ".csv"))
            utils::write.csv(miss_pentads, log_name, row.names = FALSE)
        }

        # Remove visit data from missing pentads
        visit_data <- visit_data %>%
            dplyr::anti_join(miss_pentads, by = c("Pentad", "year"))

    }

    return(list(site = site_data,
                visit = visit_data))

}
