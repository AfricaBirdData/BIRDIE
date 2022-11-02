#' Prepare Google Earth Engine visit data
#'
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambOccuR}
#' @param monitor Logical. If TRUE (default) monitoring printed messages produced
#' by `rgee` will displayed. If FALSE, only high-level messages will be displayed.
#'
#' @return
#' @export
#'
#' @examples
prepGEEVisitData <- function(config, monitor){

    # Initialize Earth Engine
    rgee::ee_check()
    rgee::ee_Initialize(drive = TRUE)

    years <- config$years

    # Load pentads to GEE -----------------------------------------------------

    # Load ABAP pentads
    pentads_sa <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa") # HARD CODED

    for(i in seq_along(years)){

        # Download SABAP data for any species -------------------------------------

        year <- years[i]

        visit <- ABAP::getAbapData(.spp_code = 6,
                                   .region_type = "country",
                                   .region = "South Africa",
                                   .years = year)


        # Annotate data with NDVI values ------------------------------------------

        # We will use the NDVI values closer to the date of the visit

        # Make spatial object and select relevant columns
        visit <- visit %>%
            dplyr::left_join(pentads_sa,
                      by = c("Pentad" = "Name")) %>%
            sf::st_sf() %>%
            dplyr::filter(!sf::st_is_empty(.)) %>%     # Remove rows without geometry
            dplyr::mutate(Date = as.character(StartDate)) %>%   # GEE doesn't like dates
            dplyr::select(CardNo, StartDate, Date, Pentad, TotalHours)

        # Upload to GEE
        ee_visit <- visit %>%
            dplyr::select(-c(StartDate, TotalHours)) %>%
            rgee::sf_as_ee(via = "getInfo")

        # Annotate with GEE NDVI
        visit_new <- ABDtools::addVarEEclosestImage(ee_feats = ee_visit,
                                                    collection = "MODIS/006/MOD13A2",
                                                    reducer = "mean",                          # We only need spatial reducer
                                                    maxdiff = 15,                              # This is the maximum time difference that GEE checks
                                                    bands = c("NDVI"),
                                                    monitor = monitor)

        visit <- visit %>%
            sf::st_drop_geometry() %>%
            dplyr::left_join(visit_new %>%
                                 sf::st_drop_geometry() %>%
                                 dplyr::select(CardNo, NDVI_mean) %>%
                                 dplyr::rename(ndvi = NDVI_mean),
                             by = c("CardNo"))

        # Update
        if(i != 1){
            visit <- rbind(visitdata, visit)
        }

        visitdata <- visit

    }


    # Prepare variables for fitting -------------------------------------------

    visitdata <- visitdata %>%
        dplyr::mutate(Date = lubridate::date(Date),
                      year = lubridate::year(Date),
                      month = lubridate::month(Date),
                      ndvi = ndvi/1e4)

    outfile <- file.path(config$out_dir, paste0("visit_dat_sa_gee_", config$years_ch, ".csv"))

    utils::write.csv(visitdata, outfile, row.names = FALSE)

    message(paste("Visits data with GEE covts saved at", outfile))

    return(visitdata)

}
