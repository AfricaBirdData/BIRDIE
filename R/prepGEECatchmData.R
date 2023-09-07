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
