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

    # Define function
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

    # Annotate with National Wetland Map --------------------------------------

    # Count pixels with wetland (wetland extension)
    out <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                   image = file.path(rgee::ee_get_assethome(), 'wetland_map_sa'),
                                   reducer = "count",
                                   monitor = monitor)

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
                                   monitor = monitor)

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

    # Last year in the GEE collection
    last_year <- 2021
    ee_collection <- rgee::ee$ImageCollection("JRC/GSW1_4/YearlyHistory")

    if(any(config$years > last_year)){

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


            # Annotate with water extension (count of pixels)
            out_count <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                                 image = ee_layer,
                                                 reducer = "count",
                                                 bands = band,
                                                 monitor = monitor,
                                                 unmask = FALSE)

            # Annotate with water recurrence (mean water class)
            out_mean <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                                image = ee_layer,
                                                reducer = "mean",
                                                bands = band,
                                                monitor = monitor,
                                                unmask = FALSE)

            # Rename because variables are named after reducer when annotated with
            # images and not with year, like when annotated with collections
            var_name <- paste0(band, "_", years_before)

            out_count <- out_count %>%
                dplyr::rename(!!var_name := paste0(band, "_count"))

            out_mean <- out_mean %>%
                dplyr::rename(!!var_name := paste0(band, "_mean"))

        } else {

            yr_range <- c(min(years_before), max(years_before))

            # Find mean for each pixel and year
            stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                                 dates = paste0(yr_range + c(0,1), "-01-01"),
                                                                 band = band,
                                                                 group_type = "year",
                                                                 groups = years_before,
                                                                 reducer = "mean",
                                                                 unmask = FALSE)

            # Find count for each pentad and year
            out_count <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "count", monitor = monitor)

            # Mean recurrence of pixels with water each year
            out_mean <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

        }

        # Create columns for years > last_year
        for(y in seq_along(years_after)){

            var_name <- paste0(band, "_", years_after[y])

            out_count <- out_count %>%
                dplyr::mutate(!!var_name := paste0(band, "_", last_year))

            out_mean <- out_mean %>%
                dplyr::mutate(!!var_name := paste0(band, "_", last_year))
        }

    } else {

        # Find mean human density for each pixel and year
        stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                             dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                             band = band,
                                                             group_type = "year",
                                                             groups = config$years,
                                                             reducer = "mean",
                                                             unmask = FALSE)

        # Find count for each pentad and year
        out_count <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "count", monitor = monitor)
        out_mean <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

    }

    # Fix covariates (it is important not to use "_" in names other than to separate the year
    out <- out %>%
        dplyr::select(Name, dplyr::starts_with("waterClass")) %>%
        dplyr::rename_with(~gsub("waterClass", "watext", .x), .cols = dplyr::starts_with("waterClass")) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")

    # Recurrence of pixels with water each year
    out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

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
    stackCollection <- ABDtools::EEcollectionToMultiband(collection = "MODIS/006/MOD13A2",
                                                         dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                         band = band,
                                                         group_type = "year",
                                                         groups = config$years,
                                                         reducer = "mean",
                                                         unmask = FALSE)

    # Find mean (mean) NDVI for each pentad and year
    out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

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

    band <- "population"

    if(any(config$years > 2020)){

        message("Using 2020 population density for years after 2020")

        years_after <- config$years[config$years > 2020]
        years_before <- config$years[!config$years > 2020]

        # Add year 2020 to retrieve info from the last year available
        if(length(years_before) == 0){
            years_before <- 2020
        }

        if(length(years_before) == 1){

            dates <- paste0(range(years_before) + c(0,1), "-01-01")

            ee_layer <- ee_collection$
                select(band)$
                filterDate(dates[1], dates[2])$
                first()

            out <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                           image = ee_layer,
                                           reducer = "mean",
                                           bands = band,
                                           monitor = monitor)

            # Rename because variables are named after reducer when annotated with
            # images and not with year, like when annotated with collections
            var_name <- paste0(band, "_", years_before)

            out <- out %>%
                dplyr::rename(!!var_name := population_mean)

        } else {

            yr_range <- c(min(years_before), max(years_before))

            # Find mean human density for each pixel and year
            stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                                 dates = paste0(yr_range + c(0,1), "-01-01"),
                                                                 band = band,
                                                                 group_type = "year",
                                                                 groups = years_before,
                                                                 reducer = "mean",
                                                                 unmask = FALSE)

            # Find mean (mean) human density for each pentad and year
            out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

        }

        # Create columns for years > 2020
        for(y in seq_along(years_after)){

            var_name <- paste0(band, "_", years_after[y])

            out <- out %>%
                dplyr::mutate(!!var_name := population_2020)
        }

    } else {

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
