#' Prepare Google Earth Engine site data
#'
#' @param config A list with pipeline configuration parameters.
#' See \code{\link{configPipeline}}
#' @param pentads An sf spatial object with the set of pentads that we want to
#' annotate with environmental covariates. See \code{\link[ABAP]{getRegionPentads}}
#' @param asset_id Character string with the name given to the object created
#' in Google Earth Engine (asset) that contains the sites in `pentads`.
#' @param upload_asset If TRUE (default), the object `pentads` will be uploaded to Google
#' Earth Engine and an asset under the name of `asset_id` will be created. If FALSE,
#' it will be assumed that an asset named after `asset_id` is already
#' present in GEE and `pentads` will not be uploaded.
#' @param monitor Logical. If TRUE (default) monitoring printed messages produced
#' by `rgee` will displayed. If FALSE, only high-level messages will be displayed.
#' @details
#' Note that some GEE layers don't have information past a certain date. At the time
#' of writing surface water layers only have information up until 2021 and human
#' population density up until 2020. We have set up the code in such a way that
#' data past the last date of the layer get annotated with the latest available
#' information. Take this into consideration for the analyses. Code should be updated
#' as more information becomes available.
#' @return
#' @export
#'
#' @examples
prepGEESiteData <- function(config, pentads, asset_id,
                            upload_asset = TRUE, monitor = TRUE){

    # Load pentads to/from to GEE ---------------------------------------------

    # Set a name for the asset
    eeid <- file.path(rgee::ee_get_assethome(), asset_id)

    # Upload to EE (if not done already)
    if(upload_asset){
        pentads %>%
            ABDtools::uploadFeaturesToEE(asset_id = eeid,
                                         load = FALSE,
                                         monitor = monitor)
    }


    # Load pentads from GEE
    ee_pentads <- rgee::ee$FeatureCollection(eeid)


    # Annotate with TerraClimate ----------------------------------------------

    message("Annotating ABAP site data with TerraClimate")

    # Define bands
    bands <- c("pr", "tmmn", "tmmx")

    # Define multi-year function
    if(length(config$years) > 1){

        f <- function(band, years, .ee_pentads, .config, .monitor){

            stackCollection <- ABDtools::EEcollectionToMultiband(collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                                 dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                                 band = band,
                                                                 group_type = "year",
                                                                 groups = years,
                                                                 reducer = "mean",
                                                                 unmask = FALSE)

            # Extract values from all bands of the image
            out <- ABDtools::addVarEEimage(.ee_pentads, stackCollection, "mean", monitor = .monitor)

            return(out)

        }

    } else {

        f <- function(band, years, .ee_pentads, .config, .monitor){

            out <- ABDtools::addVarEEcollection(ee_feats = .ee_pentads,
                                                collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                temp_reducer = "mean",
                                                spt_reducer = "mean",
                                                bands = band,
                                                unmask = FALSE,
                                                monitor = .monitor,
                                                transfer = list(via = "drive", container = "rgee_backup"))

            return(out)

        }

    }

    pentad_vars <- purrr::map(bands, ~f(.x, config$years, ee_pentads, config, monitor))

    sitedata <- Reduce("cbind", pentad_vars)

    # Fix covariates (it is important not to use "_" in names other than to separate the year
    sitedata <- sitedata %>%
        dplyr::select(!contains(".")) %>%
        dplyr::select(Name, dplyr::everything()) %>%
        dplyr::select(-id) %>%
        dplyr::rename_with(.fn =  ~gsub("pr_", "prcp_", .x),
                           .cols = dplyr::starts_with("pr_")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("tmmn"),
                                    .fns = ~.x/10),
                      dplyr::across(.cols = dplyr::starts_with("tmmx"),
                                    .fns = ~.x/10)) %>%
        sf::st_drop_geometry()

    if(length(config$years) == 1){
        sitedata <- sitedata %>%
            dplyr::rename_with(.fn =  ~gsub("_mean", paste0("_", config$years), .x),
                               .cols = dplyr::ends_with("_mean"))
    }



    # Annotate with National Wetland Map --------------------------------------

    # Count pixels with wetland (wetland extension)
    out <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                   image = file.path(rgee::ee_get_assethome(), 'wetland_map_sa'),
                                   reducer = "count",
                                   monitor = monitor,
                                   transfer = list(via = "drive", container = "rgee_backup"))

    out <- out %>%
        dplyr::rename(wetext_2018 = count) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out %>%
                             dplyr::select(Name, wetext_2018),
                         by = "Name")

    # Average wetland condition (wetland extension)
    out <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                   image = file.path(rgee::ee_get_assethome(), 'wetland_map_sa'),
                                   reducer = "mean",
                                   monitor = monitor,
                                   transfer = list(via = "drive", container = "rgee_backup"))

    sitedata <- sitedata %>%
        dplyr::left_join(out %>%
                             dplyr::rename(wetcon_2018 = mean) %>%
                             dplyr::select(Name, wetcon_2018) %>%
                             sf::st_drop_geometry(),
                         by = "Name")

    # Sites with no wetland condition are set to 0
    sitedata$wetcon_2018[is.na(sitedata$wetcon_2018)] <- 0


    # Annotate with yearly surface water occurrence --------------------------------

    message("Annotating ABAP site data with JRC surface water")

    # Number of pixels with water each year
    band <- "waterClass"

    # The JRC dataset covers up until 2021, so we annotate any years beyond this
    # with 2021 information.
    last_year <- 2021             # HARD CODED
    ee_collection <- rgee::ee$ImageCollection("JRC/GSW1_4/YearlyHistory")

    # First, annotate with water pixel count (water extension). Later with pixel
    # mean (mean water recurrence)
    reducer <- "count"

    if(any(config$years > last_year)){

        out <- addSiteVarEETimeLimit(ee_pentads, ee_collection, band, last_year,
                                     reducer=reducer, unmask=FALSE,
                                     monitor= monitor, config=config)

    } else {

        if(length(config$years) > 1){

            # Find mean human density for each pixel and year
            stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                                 dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                                 band = band,
                                                                 group_type = "year",
                                                                 groups = config$years,
                                                                 reducer = "mean",
                                                                 unmask = FALSE)

            # Find count for each pentad and year
            out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "count", monitor = monitor,
                                           transfer = list(via = "drive", container = "rgee_backup"))

        } else {

            out <- ABDtools::addVarEEcollection(ee_feats = ee_pentads,
                                                collection = ee_collection,
                                                dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                temp_reducer = "mean",
                                                spt_reducer = "count",
                                                bands = band,
                                                unmask = FALSE,
                                                monitor = .monitor,
                                                transfer = list(via = "drive", container = "rgee_backup"))

        }

    }

    # Fix covariates (it is important not to use "_" in names other than to separate the year
    out <- out %>%
        dplyr::select(Name, dplyr::starts_with("waterClass")) %>%
        dplyr::rename_with(~gsub("waterClass", "watext", .x), .cols = dplyr::starts_with("waterClass")) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")

    # We now annotate with water recurrence pixel mean (mean water recurrence)
    reducer <- "mean"

    if(any(config$years > last_year)){

        out <- addSiteVarEETimeLimit(ee_pentads, ee_collection, band, last_year,
                                     reducer=reducer, unmask=FALSE,
                                     monitor= monitor, config=config)

    } else {

        if(length(config$years) > 1){

            # Find mean human density for each pixel and year
            stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                                 dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                                 band = band,
                                                                 group_type = "year",
                                                                 groups = config$years,
                                                                 reducer = "mean",
                                                                 unmask = FALSE)

            # Find mean for each pentad and year
            out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

        } else {

            out <- ABDtools::addVarEEcollection(ee_feats = ee_pentads,
                                                collection = ee_collection,
                                                dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                temp_reducer = "mean",
                                                spt_reducer = "mean",
                                                bands = band,
                                                unmask = FALSE,
                                                monitor = .monitor,
                                                transfer = list(via = "drive", container = "rgee_backup"))

        }

    }

    # Fix covariates (it is important not to use "_" in names other than to separate the year
    out <- out %>%
        dplyr::select(Name, dplyr::starts_with("waterClass")) %>%
        dplyr::rename_with(~gsub("waterClass", "watrec", .x), .cols = dplyr::starts_with("waterClass")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("wat"),
                                    .fns = ~tidyr::replace_na(.x, 0))) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")


    # Annotate with overall surface water occurrence --------------------------

    out <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                   image = "JRC/GSW1_4/GlobalSurfaceWater",
                                   reducer = "mean",
                                   bands = "occurrence",
                                   unmask = TRUE,
                                   monitor = monitor)

    # Fix covariates (it is important not to use "_" in names other than to separate the year
    out <- out %>%
        dplyr::select(Name, watocc_ever = occurrence_mean) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")


    # Annotate with yearly NDVI -----------------------------------------------

    message("Annotating ABAP site data with MODIS NDVI")

    # Find mean NDVI for each pixel and year
    band <- "NDVI"

    if(length(config$years) > 1){

        stackCollection <- ABDtools::EEcollectionToMultiband(collection = "MODIS/006/MOD13A2",
                                                             dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                             band = band,
                                                             group_type = "year",
                                                             groups = config$years,
                                                             reducer = "mean",
                                                             unmask = FALSE)

        # Find mean (mean) NDVI for each pentad and year
        out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

    } else {

        out <- ABDtools::addVarEEcollection(ee_feats = ee_pentads,
                                            collection = "MODIS/006/MOD13A2",
                                            dates = paste0(config$year_range + c(0,1), "-01-01"),
                                            temp_reducer = "mean",
                                            spt_reducer = "mean",
                                            bands = band,
                                            unmask = FALSE,
                                            monitor = monitor,
                                            transfer = list(via = "drive", container = "rgee_backup"))

    }

    # Fix covariates (it is important not to use "_" in names other than to separate the year
    out <- out %>%
        dplyr::select(Name, dplyr::starts_with("NDVI")) %>%
        dplyr::rename_with(~gsub("NDVI", "ndvi", .x), .cols = dplyr::starts_with("NDVI")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("ndvi"),
                                    .fns = ~.x/10000)) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")


    # Annotate with yearly human population density ---------------------------

    message("Annotating ABAP site data with WorldPop")

    # We need to subset the image collection to South Africa
    ee_collection <- rgee::ee$ImageCollection("WorldPop/GP/100m/pop")
    ee_collection <- ee_collection$filterMetadata('country','equals','ZAF')

    # We will use the last year present in the data set
    last_year <- 2020              # HARD CODED
    band <- "population"
    reducer <- "mean"

    if(any(config$years > last_year)){

        out <- addSiteVarEETimeLimit(ee_pentads, ee_collection, band, last_year,
                                     reducer=reducer, unmask=FALSE,
                                     monitor= monitor, config=config)

    } else {

        if(length(config$years) > 1){

            # Find mean human density for each pixel and year
            stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                                 dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                                 band = band,
                                                                 group_type = "year",
                                                                 groups = config$years,
                                                                 reducer = "mean",
                                                                 unmask = FALSE)

            # Find mean (mean) human density for each pentad and year
            out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

        } else {

            out <- ABDtools::addVarEEcollection(ee_feats = ee_pentads,
                                                collection = ee_collection,
                                                dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                temp_reducer = "mean",
                                                spt_reducer = "mean",
                                                bands = band,
                                                unmask = FALSE,
                                                monitor = monitor,
                                                transfer = list(via = "drive", container = "rgee_backup"))

        }

    }

    # Fix covariates (it is important not to use "_" in names other than to separate the year
    out <- out %>%
        dplyr::select(Name, dplyr::starts_with("population")) %>%
        dplyr::rename_with(~gsub("population", "hum.km2", .x), .cols = dplyr::starts_with("population")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("hum.km2"),
                                    .fns = ~.x*100)) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")


    # Add geometries -----------------------------------------------------------

    sitedata <- sitedata %>%
        dplyr::left_join(dplyr::select(pentads, Name), by = "Name") %>%
        sf::st_sf()

    # Add centroids
    # sf::st_agr(sitedata) = "constant"
    cc <- sf::st_centroid(sitedata) %>%
        sf::st_coordinates() %>%
        as.data.frame()

    sitedata <- sitedata %>%
        dplyr::mutate(lon = cc$X,
                      lat = cc$Y)


    # Annotate with distance to coast -----------------------------------------

    message("Annotating ABAP site data with distance to coast and elevation")

    # Create bounding box
    sabox <- sf::st_bbox(sitedata) + c(-1, -1, 3, 1)

    # Download coastline for the World and crop
    coast <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf")
    #sf::st_agr(coast) = "constant"
    coast <- coast %>%
        sf::st_crop(sabox)

    # Find distances in kilometers
    sitedata$dist_coast <- as.numeric(sf::st_distance(sitedata, coast[1,]))/1000 # in kilometers


    # Annotate with elevation -------------------------------------------------

    out <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                   image = "MERIT/DEM/v1_0_3",
                                   reducer = "mean",
                                   bands = "dem",
                                   unmask = FALSE,
                                   monitor = monitor)

    # Fix covariates (it is important not to use "_" in names other than to separate the year
    out <- out %>%
        dplyr::select(Name, elev = dem_mean) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        sf::st_drop_geometry() %>%
        dplyr::left_join(out, by = "Name")

    return(sitedata)

}


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
                    dplyr::distinct(Date, Pentad, .keep_all = TRUE) %>%
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


#' Prepare Google Earth Engine data for catchments
#'
#' @param sp_code SAFRING reference number of the species we want to analyze.
#' @param catchment An sf object with the polygons defining the catchments
#' to be annotated.
#' @param config A list with pipeline configuration parameters.
#' See \code{\link{configPipeline}}
#' @param monitor Logical. If TRUE (default) monitoring printed messages produced
#' by `rgee` will displayed. If FALSE, only high-level messages will be displayed.
#'
#' @details It is assumed that there is an asset on `rgee::ee_get_assethome()`
#' named 'quin_catchm' that has the polygons defining the quinary catchments.
#' Note that we add an extra year before the start of the series. This is because
#' summer waterbird populations should be affected by the conditions in the previous
#' year rather than by conditions in the following year.
#' @return
#' @export
#'
#' @examples
prepGEECatchmData <- function(sp_code, catchment, config, monitor = TRUE){

    # Note that we add an extra year before the start of the series. This is because
    # summer waterbird populations should be affected by the conditions in the previous
    # year rather than by conditions in the following year.
    config$years <- c(min(config$years) - 1, config$years)
    config$year_range <- range(config$years)


    # Load quinaries to GEE -----------------------------------------------------

    # Set a name for the asset
    eeCatchm_id <- file.path(rgee::ee_get_assethome(), 'quin_catchm')

    # Load catchm from GEE
    ee_catchm <- rgee::ee$FeatureCollection(eeCatchm_id)


    # Annotate with TerraClimate ----------------------------------------------

    message("Annotating catchment data with TerraClimate")

    # Define bands
    bands <- c("pr", "tmmn", "tmmx", "pdsi")

    # Define function
    f <- function(band, years, .ee_catchm, .config, .monitor){

        stackCollection <- ABDtools::EEcollectionToMultiband(collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                             dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                             band = band,
                                                             group_type = "year",
                                                             groups = years,
                                                             reducer = "mean",
                                                             unmask = FALSE)

        # Extract values from all bands of the image
        out <- ABDtools::addVarEEimage(.ee_catchm, stackCollection, "mean", monitor = .monitor)

        return(out)

    }

    catchm_vars <- purrr::map(bands, ~f(.x, config$years, ee_catchm, config, monitor))

    sitedata <- Reduce("cbind", catchm_vars)

    # Fix covariates
    sitedata <- sitedata %>%
        dplyr::select(!contains(".")) %>%
        dplyr::select(UNIT_ID, dplyr::everything()) %>%
        dplyr::select(-id) %>%
        dplyr::rename_with(.fn =  ~gsub("pr_", "prcp_", .x),
                           .cols = dplyr::starts_with("pr_")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("tmmn"),
                                    .fns = ~.x/10),
                      dplyr::across(.cols = dplyr::starts_with("tmmx"),
                                    .fns = ~.x/10)) %>%
        sf::st_drop_geometry()


    # Annotate with National Wetland Map --------------------------------------

    message("Annotating catchment data with National Wetland Map")

    # Count pixels with wetland (wetland extension)
    out <- ABDtools::addVarEEimage(ee_feats = ee_catchm,
                                   image = file.path(rgee::ee_get_assethome(), 'wetland_map_sa'),
                                   reducer = "count",
                                   monitor = TRUE)

    out <- out %>%
        dplyr::mutate(area.km2 = as.numeric(sf::st_area(.) %>% units::set_units(km^2)),
                      wetden = count/area.km2,
                      log_wetden = log(wetden+1)) %>%
        dplyr::rename(wetext = count) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out %>%
                             dplyr::select(UNIT_ID, wetext, wetden, area.km2),
                         by = "UNIT_ID")

    # Average wetland condition (wetland extension)
    out <- ABDtools::addVarEEimage(ee_feats = ee_catchm,
                                   image = file.path(rgee::ee_get_assethome(), 'wetland_map_sa'),
                                   reducer = "mean",
                                   monitor = TRUE)

    sitedata <- sitedata %>%
        dplyr::left_join(out %>%
                             dplyr::rename(wetcon = mean) %>%
                             dplyr::select(UNIT_ID, wetcon) %>%
                             sf::st_drop_geometry(),
                         by = "UNIT_ID")

    # Sites with no wetland condition are set to 0
    sitedata$wetcon[is.na(sitedata$wetcon)] <- 0

    # Annotate with yearly surface water occurrence --------------------------------

    message("Annotating catchment data with JRC surface water")

    # This data set runs up to 2021 so we can't annotate anything further
    jrc_years <- config$years[config$years < 2022]
    jrc_year_range <- range(jrc_years)
    years_miss <- config$years[config$years > 2021]

    # Number of pixels with water each year
    band <- "waterClass"
    stackCollection <- ABDtools::EEcollectionToMultiband(collection = "JRC/GSW1_4/YearlyHistory",
                                                         dates = paste0(jrc_year_range + c(0,1), "-01-01"),
                                                         band = band,
                                                         group_type = "year",
                                                         groups = jrc_years,
                                                         reducer = "mean",
                                                         unmask = FALSE)

    out <- ABDtools::addVarEEimage(ee_catchm, stackCollection, "count", monitor = monitor)

    # Fix covariates
    out <- out %>%
        dplyr::select(UNIT_ID, dplyr::starts_with("waterClass")) %>%
        dplyr::rename_with(~gsub("waterClass", "watext", .x), .cols = dplyr::starts_with("waterClass")) %>%
        sf::st_drop_geometry()

    # Add NAs to those years that have no data
    for(i in seq_along(years_miss)){
        vv <- paste0("watext_", years_miss[i])
        out <- out %>%
            dplyr::mutate(!!vv := NA)
    }

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "UNIT_ID")

    # Recurrence of pixels with water each year
    out <- ABDtools::addVarEEimage(ee_catchm, stackCollection, "mean", monitor = monitor)

    # Fix covariates
    out <- out %>%
        dplyr::select(UNIT_ID, dplyr::starts_with("waterClass")) %>%
        dplyr::rename_with(~gsub("waterClass", "watrec", .x), .cols = dplyr::starts_with("waterClass")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("wat"),
                                    .fns = ~tidyr::replace_na(.x, 0))) %>%
        sf::st_drop_geometry()

    # Add NAs to those years that have no data
    for(i in seq_along(years_miss)){
        vv <- paste0("watrec_", years_miss[i])
        out <- out %>%
            dplyr::mutate(!!vv := NA)
    }

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "UNIT_ID")


    # Annotate with overall surface water occurrence --------------------------

    out <- ABDtools::addVarEEimage(ee_feats = ee_catchm,
                                   image = "JRC/GSW1_4/GlobalSurfaceWater",
                                   reducer = "mean",
                                   bands = "occurrence",
                                   unmask = TRUE,
                                   monitor = monitor)

    # Fix covariates
    out <- out %>%
        dplyr::select(UNIT_ID, watocc_ever = occurrence_mean) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "UNIT_ID")


    # Annotate with yearly NDVI -----------------------------------------------

    message("Annotating catchment data with MODIS NDVI")

    # Find mean water presence for each pixel and year
    band <- "NDVI"

    if(any(config$years < 2000)){

        years_before <- config$years[config$years < 2000]
        years_after <- config$years[!config$years < 2000]

        stackCollection <- ABDtools::EEcollectionToMultiband(collection = "NASA/GIMMS/3GV0",
                                                             dates = paste0(range(years_before) + c(0,1), "-01-01"),
                                                             band = tolower(band),
                                                             group_type = "year",
                                                             groups = years_before,
                                                             reducer = "mean",
                                                             unmask = FALSE)

        # Find mean (mean) NDVI for each catchm and year
        out1 <- ABDtools::addVarEEimage(ee_catchm, stackCollection, "mean", monitor = monitor)

        # Prepare for joining
        out1 <- out1 %>%
            sf::st_drop_geometry() %>%
            dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("ndvi"),
                                        .fns = ~.x*10000))

    } else {
        out1 <- data.frame()
        years_after <- config$years
    }

    stackCollection <- ABDtools::EEcollectionToMultiband(collection = "MODIS/006/MOD13A2",
                                                         dates = paste0(range(years_after) + c(0,1), "-01-01"),
                                                         band = band,
                                                         group_type = "year",
                                                         groups = years_after,
                                                         reducer = "mean",
                                                         unmask = FALSE)

    # Find mean (mean) NDVI for each catchm and year
    out <- ABDtools::addVarEEimage(ee_catchm, stackCollection, "mean", monitor = monitor)

    # Fix covariates
    out <- out1 %>%
        dplyr::left_join(out %>%
                             sf::st_drop_geometry() %>%
                             dplyr::select(UNIT_ID, dplyr::starts_with("NDVI")),
                         by = "UNIT_ID") %>%
        dplyr::rename_with(~gsub("NDVI", "ndvi", .x), .cols = dplyr::starts_with("NDVI")) %>%
        dplyr::select(UNIT_ID, dplyr::starts_with("NDVI")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("ndvi"),
                                    .fns = ~.x/10000))

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "UNIT_ID")


    # Add geometries -----------------------------------------------------------

    sitedata <- sitedata %>%
        dplyr::left_join(dplyr::select(catchment, UNIT_ID), by = "UNIT_ID") %>%
        sf::st_sf()

    # Add centroids
    # sf::st_agr(sitedata) = "constant"
    cc <- sf::st_centroid(sitedata) %>%
        sf::st_coordinates() %>%
        as.data.frame()

    sitedata <- sitedata %>%
        dplyr::mutate(lon = cc$X,
                      lat = cc$Y)


    # Annotate with distance to coast -----------------------------------------

    message("Annotating catchment data with distance to coast and elevation")

    # Create bounding box
    sabox <- sf::st_bbox(sitedata) + c(-1, -1, 3, 1)

    # Download coastline for the World and crop
    coast <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf")
    #sf::st_agr(coast) = "constant"
    coast <- coast %>%
        sf::st_crop(sabox)

    # Find distances in kilometers
    sitedata$dist_coast <- as.numeric(sf::st_distance(sitedata, coast[1,]))/1000 # in kilometers


    # Annotate with elevation -------------------------------------------------

    out <- ABDtools::addVarEEimage(ee_feats = ee_catchm,
                                   image = "MERIT/DEM/v1_0_3",
                                   reducer = "mean",
                                   bands = "dem",
                                   unmask = FALSE,
                                   monitor = monitor)

    # Fix covariates
    out <- out %>%
        dplyr::select(UNIT_ID, elev = dem_mean) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        sf::st_drop_geometry() %>%
        dplyr::left_join(out, by = "UNIT_ID")


    # Round results
    sitedata <- sitedata %>%
        dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), round, 3))


    # Return ------------------------------------------------------------------

    return(sitedata)

}


#' Annotate site data with Google Earth Engine with a time limit
#'
#' @param ee_feats A feature collection with the sites we want to annotate.
#' We should have uploaded an sf object with the sites to GEE, previously.
#' See \code{\link{uploadFeaturesToEE}}
#' @param ee_collection A GEE collection produced with \code{ee$ImageCollection()}.
#' See \href{https://developers.google.com/earth-engine/datasets/catalog}{GEE catalog}.
#' @param band The band of the collection we are going to use to annotate our data
#' @param last_year The last year provided by the GEE image collection
#' @param reducer The reducer we will use to summarize the images of the image collection
#' @param unmask GEE masks missing values, which means they are not used for
#' computing means, counts, etc. Sometimes we might want to avoid this behaviour
#' and use 0 instead of NA. If so, set unmask to TRUE.
#' @param monitor Logical. If TRUE (default) monitoring messages produced
#' by `rgee` will displayed. If FALSE, only high-level messages will be displayed.
#' @param config A list with pipeline configuration parameters.
#' See \code{\link{configPipeline}}
#' @details
#' The function annotates data with the corresponding year if it is present in the
#' GEE image collection. Years after the last year present in 'config$years' are annotated
#' with the last year present in the image collection.
#' @return
#' @export
#'
#' @examples
addSiteVarEETimeLimit <- function(ee_feats, ee_collection, band, last_year, reducer,
                                  unmask, monitor, config){

    message(paste("Using", last_year, "for", band, "for years after", last_year))

    years_after <- config$years[config$years > last_year]
    years_before <- config$years[!config$years > last_year]

    # Add year "last year" to retrieve info from the last year available
    if(length(years_before) == 0){
        years_before <- last_year
    }

    if(length(years_before) == 1){

        dates <- paste0(range(years_before) + c(0,1), "-01-01")

        ee_layer <- ee_collection$
            select(band)$
            filterDate(dates[1], dates[2])$
            first()


        # Annotate with image
        out <- ABDtools::addVarEEimage(ee_feats = ee_feats,
                                       image = ee_layer,
                                       reducer = reducer,
                                       bands = band,
                                       monitor = monitor,
                                       unmask = unmask,
                                       transfer = list(via = "drive", container = "rgee_backup"))

        # Rename because variables are named after reducer when annotated with
        # images and not with year, like when annotated with collections
        var_name <- paste0(band, "_", years_before)

        out <- out %>%
            dplyr::rename(!!var_name := .data[[paste0(band, "_", reducer)]] )

    } else {

        yr_range <- c(min(years_before), max(years_before))

        # Make multiband image
        stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                             dates = paste0(yr_range + c(0,1), "-01-01"),
                                                             band = band,
                                                             group_type = "year",
                                                             groups = years_before,
                                                             reducer = "mean",     # HARD CODED. Temporal groups always get a mean. Check that this makes sense
                                                             unmask = unmask)

        # Annotate with image
        out <- ABDtools::addVarEEimage(ee_feats, stackCollection, reducer, monitor = monitor,
                                       transfer = list(via = "drive", container = "rgee_backup"))

    }

    # Create columns for years > last_year
    for(y in seq_along(years_after)){

        var_name <- paste0(band, "_", years_after[y])

        out <- out %>%
            dplyr::mutate(!!var_name := .data[[paste0(band, "_", last_year)]])

    }

    return(out)

}

