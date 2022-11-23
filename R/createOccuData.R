#' Create data for fitting an occupancy model
#'
#' @description This function takes generic occupancy visit and site data
#' without species observations, adds detection data for a species in a given
#' year, runs some checks and formats the output.
#' @param sp_code SAFRING_No of the species of interest as extracted from ABAP.
#' Ignored if download set to FALSE.
#' @param years A numeric vector with elements corresponding to the years we
#' want data for. Ignored if download set to FALSE.
#' @param site_data A dataframe with occupancy site data.
#' @param visit_data A dataframe with occupancy visit data.
#' @param force_abap_dwld Indicates whether ABAP data must be downloaded for the species
#' and years indicated by 'sp_code' and 'years'. If TRUE, data will be downloaded
#' from ABAP once per session and cached in a temp file. After this the cached
#' file will be used, unless download is set to FALSE, in which case data will
#' be downloaded regardless of the cached file.
#' @param config A list of configuration parameters see \link{configPreambOccu}
#'
#' @return A list containing two data frames: one with site data and one with
#' visit data.
#' @export
#'
#' @examples
createOccuData <- function(sp_code, years,
                           site_data, visit_data,
                           force_abap_dwld = FALSE, config){


    # Cache file name
    cachefile <- file.path(tempdir(), paste(c(sp_code, years, ".rds"), collapse = "_"))

    if(!file.exists(cachefile) | force_abap_dwld){

        # Download species detection
        message("Downloading from SABAP")

        sp_detect <- ABAP::getAbapData(.spp_code = sp_code,
                                       .region_type = "country",
                                       .region = "South Africa",
                                       .years = years)

        # Save to cache
        saveRDS(sp_detect, cachefile)

    } else {

        message("Using cached file")
        sp_detect <- readRDS(cachefile)

    }

    # Other variables and remove those rows with NA values
    site_data <- site_data %>%
        dplyr::select(Pentad = Name, lon, lat, watocc_ever, dist_coast, elev, dplyr::ends_with(match = as.character(config$years))) %>%
        tidyr::drop_na()   # I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE

    # Subset visits that correspond to the years of interest
    visit_data <- visit_data %>%
        dplyr::filter(year %in% config$years)

    # Add detections to visit data
    visit_data <- visit_data %>%
        dplyr::mutate(StartDate = as.Date(StartDate)) %>%
        dplyr::left_join(sp_detect %>%
                             dplyr::select(CardNo, StartDate, Pentad, obs = Spp) %>%
                             dplyr::mutate(obs = ifelse(obs == "-", 0, 1)),
                         by = c("CardNo", "StartDate", "Pentad"))


    if(any(is.na(visit_data$obs))){
        warning("NA found in detection data")
        visit_data <- visit_data %>%
            dplyr::filter(!is.na(obs))
    }


    # Format site data --------------------------------------------------------

    # First, long format for site variables and years
    # HARD CODED! Check that these are the variables that don't change over time
    fixed_vars <- c("Pentad", "lon", "lat", "site", "watocc_ever", "dist_coast", "elev")

    site_data <- site_data %>%
        BIRDIE::gatherYearFromVars(vars = setdiff(names(.), fixed_vars), sep = "_")


    # Create other covariates -------------------------------------------------

    message("Creating additional covariates...")

    # Load South African pentads
    sa_pentads <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa")

    # Create a detection covariate for CWAC sites
    cwac_sites <- CWAC::listCwacSites(.region_type = "country", .region = "South Africa") %>%
        sf::st_as_sf(coords = c("X", "Y"), crs = 4326, remove = FALSE)

    # Make sites a spatial object
    site_data <- site_data %>%
        dplyr::mutate(Name = Pentad) %>%
        dplyr::left_join(sa_pentads %>%
                             dplyr::select(Name),
                         by = "Name") %>%
        sf::st_sf()

    cwac_intsc <- sf::st_intersects(site_data, cwac_sites)
    site_data <- site_data %>%
        dplyr::mutate(cwac = sapply(cwac_intsc, length) != 0) %>%
        sf::st_drop_geometry()

    # Create other variables
    site_data <- site_data %>%
        dplyr::mutate(tdiff = tmmx - tmmn,
                      log_dist_coast = log(dist_coast + 1),
                      log_watext = log(watext + 1),
                      log_watrec = log(watrec + 1))


    # Save to disc ------------------------------------------------------------
    return(list(site = site_data,
                visit = visit_data))

}
