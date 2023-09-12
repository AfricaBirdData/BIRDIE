#' Prepare Google Earth Engine visit data
#'
#' @param config A list with pipeline configuration parameters.
#' See \code{\link{configPipeline}}
#' @param visits An sf spatial object with the visit data that we want to
#' annotate with environmental covariates. Note that we need these data to be
#' in a spatial format, so we should probably join them with pentad data.
#' See \code{\link[ABAP]{getRegionPentads}}
#' @param asset_id Character string with the name given to the object created
#' in Google Earth Engine (asset) that contains the sites in `visits`.
#' @param upload_asset If TRUE (default), the object `visits` will be uploaded to Google
#' Earth Engine and an asset under the name of `asset_id` will be created. If FALSE,
#' it will be assumed that an asset named after `asset_id` is already present in GEE
#' and `visits` will not be uploaded.
#' @param monitor Logical. If TRUE (default) monitoring printed messages produced
#' by `rgee` will displayed. If FALSE, only high-level messages will be displayed.
#' @details
#' Note that some GEE layers don't have information past a certain date. At the time
#' of writing surface water layers only have information up until 2021 and human
#' population density up until 2020. We have set up the code in such a way that
#' visits past the last date of the layer get annotated with the latest available
#' information. Take this into consideration for the analyses. Code should be updated
#' as more information becomes available.
#'
#' @return
#' @export
#'
#' @examples
prepGEEVisitData <- function(config, visits, asset_id,
                             upload_asset = TRUE, monitor = TRUE){

    # Load pentads to GEE -----------------------------------------------------

    # We need annotate data year by year. Otherwise object are too big to send to GEE
    years <- config$years

    for(i in seq_along(years)){

        # Download SABAP data for any species -------------------------------------

        year_sel <- years[i]

        visit <- visits %>%
            dplyr::filter(year == year_sel)

        # Save geometry for later
        gm <- visit %>%
            dplyr::distinct(Pentad, geometry)

        # Set a name for the asset
        eeid <- file.path(rgee::ee_get_assethome(), asset_id)

        # Upload to EE (if not done already)
        if(upload_asset){
            visit %>%
                dplyr::select(-c(StartDate, TotalHours)) %>%
                dplyr::distinct(Date, Pentad, .keep_all = TRUE) %>%
                ABDtools::uploadFeaturesToEE(asset_id = eeid,
                                             load = FALSE,
                                             monitor = monitor)
        }

        # Load pentads from GEE
        ee_visit <- rgee::ee$FeatureCollection(eeid)


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

        # We need to be joining by pentad and date
        visit <- visit %>%
            sf::st_drop_geometry() %>%
            dplyr::left_join(visit_ndvi %>%
                                 sf::st_drop_geometry() %>%
                                 dplyr::select(Pentad, Date, NDVI_mean) %>%
                                 dplyr::mutate(Date = as.character(Date)) %>%
                                 dplyr::rename(ndvi = NDVI_mean),
                             by = c("Pentad", "Date"))
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
                                 dplyr::select(Pentad, Date, pr_mean) %>%
                                 dplyr::mutate(Date = as.character(Date)) %>%
                                 dplyr::rename(prcp = pr_mean),
                             by = c("Pentad", "Date"))

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
                                 dplyr::select(Pentad, Date, tmmn_mean) %>%
                                 dplyr::mutate(Date = as.character(Date)) %>%
                                 dplyr::rename(tmmn = tmmn_mean),
                             by = c("Pentad", "Date"))

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
                                 dplyr::select(Pentad, Date, tmmx_mean) %>%
                                 dplyr::mutate(Date = as.character(Date)) %>%
                                 dplyr::rename(tmmx = tmmx_mean),
                             by = c("Pentad", "Date"))

        rm(visit_tmmx)

        # Annotate data with human population -------------------------------------

        message("Annotating ABAP visit data with WorldPop")

        # Upload to EE with new dates if necessary. We will use data from 2020
        # for those visits that occur later on
        if(year_sel > 2020){

            message("Using 2020 population density for years after 2020")

            # make sure dates match 2020 images. There is only one image for each year,
            # so we set Date to be the beginning of the year 2020
            visit_2020 <- visit %>%
                dplyr::mutate(Date = "2020-02-01")

            # Add geometry
            visit_2020 <- visit_2020 %>%
                dplyr::left_join(gm, by = "Pentad") %>%
                sf::st_sf()

            # Upload asset
            if(upload_asset){
                visit_2020 %>%
                    dplyr::select(year, CardNo, Date, Pentad, ObserverNo) %>%
                    dplyr::distinct(Date, Pentad) %>%
                    ABDtools::uploadFeaturesToEE(asset_id = eeid,
                                                 load = FALSE,
                                                 monitor = monitor)
            } else {
                warning("Cannot annotate visits with data newer than 2020. Setting upload asset = TRUE would help, but read documentation first.")
            }


        }

        # Load pentads from GEE
        ee_visit <- rgee::ee$FeatureCollection(eeid)

        # We need to subset the image collection to South Africa
        ee_collection <- rgee::ee$ImageCollection("WorldPop/GP/100m/pop")
        ee_collection <- ee_collection$filterMetadata('country','equals','ZAF')

        visit_pop <- ABDtools::addVarEEclosestImage(ee_feats = ee_visit,
                                                    collection = ee_collection,
                                                    reducer = "mean",                          # We only need spatial reducer
                                                    maxdiff = 200,                              # This is the maximum time difference that GEE checks
                                                    bands = c("population"),
                                                    monitor = monitor)

        visit <- visit %>%
            sf::st_drop_geometry() %>%
            dplyr::left_join(visit_pop %>%
                                 sf::st_drop_geometry() %>%
                                 dplyr::select(CardNo, population_mean) %>%
                                 dplyr::rename(hum.km2 = population_mean),
                             by = c("CardNo"))

        rm(visit_pop)

        # Update
        if(i != 1){
            visit <- rbind(visitdata, visit)
        }

        visitdata <- visit

    }


    # Prepare variables for fitting -------------------------------------------

    # Fix covariates (it is important not to use "_" in names other than to separate the year)
    visitdata <- visitdata %>%
        dplyr::mutate(Date = lubridate::date(Date),
                      year = lubridate::year(Date),
                      month = lubridate::month(Date),
                      ndvi = ndvi/1e4,
                      tmmn = tmmn/10,
                      tmmx = tmmx/10,
                      hum.km2 = hum.km2*100)

    return(visitdata)

}
