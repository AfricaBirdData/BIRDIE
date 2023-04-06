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
prepGEEVisitData <- function(config, monitor = TRUE){

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

        visit <- ABAP::getAbapData(.spp_code = 6,  # this species should not matter - visits are the same for all species
                                   .region_type = "country",
                                   .region = "South Africa",
                                   .years = year)

        # Make spatial object and select relevant columns
        visit <- visit %>%
            dplyr::left_join(pentads_sa,
                             by = c("Pentad" = "Name")) %>%
            sf::st_sf() %>%
            dplyr::filter(!sf::st_is_empty(.)) %>%     # Remove rows without geometry
            dplyr::mutate(Date = as.character(StartDate)) %>%   # GEE doesn't like dates
            dplyr::select(CardNo, StartDate, Date, Pentad, TotalHours, ObserverNo)

        # Upload to GEE
        ee_visit <- visit %>%
            dplyr::select(-c(StartDate, TotalHours)) %>%
            rgee::sf_as_ee(via = "getInfo")


        # Annotate data with NDVI values ------------------------------------------

        # We will use the NDVI values closer to the date of the visit

        # Annotate with GEE NDVI
        message("Annotating ABAP visit data with NDVI")

        visit_ndvi <- ABDtools::addVarEEclosestImage(ee_feats = ee_visit,
                                                     collection = "MODIS/006/MOD13A2",
                                                     reducer = "mean",                          # We only need spatial reducer
                                                     maxdiff = 15,                              # This is the maximum time difference that GEE checks
                                                     bands = c("NDVI"),
                                                     monitor = monitor)

        visit <- visit %>%
            sf::st_drop_geometry() %>%
            dplyr::left_join(visit_ndvi %>%
                                 sf::st_drop_geometry() %>%
                                 dplyr::select(CardNo, NDVI_mean) %>%
                                 dplyr::rename(ndvi = NDVI_mean),
                             by = c("CardNo"))
        rm(visit_ndvi)


        # Annotate with GEE TerraClimate ------------------------------------------

        # Annotate with GEE TerraClimate precipitation

        message("Annotating ABAP visit data with TerraClimate")

        visit_prcp <- ABDtools::addVarEEclosestImage(ee_feats = ee_visit,
                                                     collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                     reducer = "mean",                          # We only need spatial reducer
                                                     maxdiff = 20,                              # This is the maximum time difference that GEE checks
                                                     bands = c("pr"),
                                                     monitor = monitor)

        visit <- visit %>%
            sf::st_drop_geometry() %>%
            dplyr::left_join(visit_prcp %>%
                                 sf::st_drop_geometry() %>%
                                 dplyr::select(CardNo, pr_mean) %>%
                                 dplyr::rename(prcp = pr_mean),
                             by = c("CardNo"))

        rm(visit_prcp)

        # Annotate with GEE TerraClimate minimum temperature
        visit_tmmn <- ABDtools::addVarEEclosestImage(ee_feats = ee_visit,
                                                     collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                     reducer = "mean",                          # We only need spatial reducer
                                                     maxdiff = 20,                              # This is the maximum time difference that GEE checks
                                                     bands = c("tmmn"),
                                                     monitor = monitor)

        visit <- visit %>%
            sf::st_drop_geometry() %>%
            dplyr::left_join(visit_tmmn %>%
                                 sf::st_drop_geometry() %>%
                                 dplyr::select(CardNo, tmmn_mean) %>%
                                 dplyr::rename(tmmn = tmmn_mean),
                             by = c("CardNo"))

        rm(visit_tmmn)

        # Annotate with GEE TerraClimate maximum temperature
        visit_tmmx <- ABDtools::addVarEEclosestImage(ee_feats = ee_visit,
                                                     collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                     reducer = "mean",                          # We only need spatial reducer
                                                     maxdiff = 20,                              # This is the maximum time difference that GEE checks
                                                     bands = c("tmmx"),
                                                     monitor = monitor)

        visit <- visit %>%
            sf::st_drop_geometry() %>%
            dplyr::left_join(visit_tmmx %>%
                                 sf::st_drop_geometry() %>%
                                 dplyr::select(CardNo, tmmx_mean) %>%
                                 dplyr::rename(tmmx = tmmx_mean),
                             by = c("CardNo"))

        rm(visit_tmmx)

        # Annotate data with human population -------------------------------------

        message("Annotating ABAP visit data with WorldPop")

        visit_pop <- ABDtools::addVarEEclosestImage(ee_feats = ee_visit,
                                                    collection = "WorldPop/GP/100m/pop",
                                                    reducer = "mean",                          # We only need spatial reducer
                                                    maxdiff = 15,                              # This is the maximum time difference that GEE checks
                                                    bands = c("population"),
                                                    monitor = monitor)

        visit <- visit %>%
            sf::st_drop_geometry() %>%
            dplyr::left_join(visit_pop %>%
                                 sf::st_drop_geometry() %>%
                                 dplyr::select(CardNo, population_mean) %>%
                                 dplyr::rename(hum_km2 = population_mean),
                             by = c("CardNo"))

        rm(visit_pop)

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
                      ndvi = ndvi/1e4,
                      tmmn = tmmn/10,
                      tmmx = tmmx/10,
                      hum_km2 = hum_km*100)


    outfile <- file.path(config$out_dir, paste0("visit_dat_sa_gee_", config$years_ch, ".csv"))

    utils::write.csv(visitdata, outfile, row.names = FALSE)

    message(paste("Visits data with GEE covts saved at", outfile))

    return(visitdata)

}
