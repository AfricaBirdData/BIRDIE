#' Prepare Google Earth Engine site data
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
prepGEESiteData <- function(config, monitor = monitor){

    # Initialize Earth Engine
    rgee::ee_check()
    rgee::ee_Initialize(drive = TRUE)


    # Load pentads to GEE -----------------------------------------------------

    # Load ABAP pentads
    pentads_sa <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa") # HARD CODED

    # Set a name for the asset
    eeid <- file.path(rgee::ee_get_assethome(), 'birdie_pentads_sa')

    # Upload to EE (if not done already)
    # ee_pentads <- pentads_sa %>%
    #     ABDtools::uploadFeaturesToEE(asset_id = eeid,
    #                                  load = TRUE,
    #                                  monitor = monitor)

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
    stackCollection <- ABDtools::EEcollectionToMultiband(collection = "JRC/GSW1_4/YearlyHistory",
                                                         dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                         band = band,
                                                         group_type = "year",
                                                         groups = config$years,
                                                         reducer = "mean",
                                                         unmask = FALSE)

    out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "count", monitor = monitor)

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

    # Find mean human density for each pixel and year
    band <- "population"
    stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                         dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                         band = band,
                                                         group_type = "year",
                                                         groups = config$years,
                                                         reducer = "mean",
                                                         unmask = FALSE)

    # Find mean (mean) human density for each pentad and year
    out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, "mean", monitor = monitor)

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
        dplyr::left_join(dplyr::select(pentads_sa, Name), by = "Name") %>%
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

    outfile <- file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv"))

    utils::write.csv(sitedata, outfile, row.names = FALSE)

    message(paste("Site data with GEE covts saved at", outfile))

    return(sitedata)

}
