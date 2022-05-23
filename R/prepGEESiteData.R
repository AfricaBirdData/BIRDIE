#' Prepare Google Earth Engine site data
#'
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambOccuR}
#'
#' @return
#' @export
#'
#' @examples
prepGEESiteData <- function(config){

    # Initialize Earth Engine
    rgee::ee_check()
    rgee::ee_Initialize(drive = TRUE)


    # Load pentads to GEE -----------------------------------------------------

    # Load ABAP pentads
    pentads_sa <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa") # HARD CODED

    # Set a name for the asset
    eeid <- file.path(rgee::ee_get_assethome(), 'pentads')

    # Upload to EE (if not done already). Delete first
    # assets <- rgee::ee_manage_assetlist(rgee::ee_get_assethome())
    #
    # if(eeid %in% assets$ID){
    #     rgee::ee_manage_delete(eeid, quiet = FALSE, strict = TRUE)
    # }

    # ABAP::uploadPentadsToEE(pentads = dplyr::select(pentads_sa, Name),
                            # asset_id = eeid,
                            # load = FALSE)

    # Load pentads from GEE
    ee_pentads <- rgee::ee$FeatureCollection(eeid)


    # Annotate with TerraClimate ----------------------------------------------

    # Define bands
    bands <- c("pr", "tmmn", "tmmx")

    # Define function
    f <- function(band, years, .ee_pentads, .config){

        stackCollection <- ABAP::EEcollectionToMultiband(collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                         dates = paste0(.config$year_range + c(0,1), "-01-01"),
                                                         band = band,
                                                         group_type = "year",
                                                         groups = years,
                                                         reducer = "mean",
                                                         unmask = FALSE)

        # Extract values from all bands of the image
        out <- ABAP::addVarEEimage(.ee_pentads, stackCollection, "mean")

        out <- out %>%
            dplyr::rename_with(~gsub("X", band, .x), .cols = dplyr::starts_with("X"))

        return(out)

    }

    pentad_vars <- purrr::map(bands, ~f(.x, config$years, ee_pentads, config))

    sitedata <- Reduce("cbind", pentad_vars)

    # Fix covariates
    sitedata <- sitedata %>%
        dplyr::select(!contains(".")) %>%
        dplyr::select(Name, dplyr::everything()) %>%
        dplyr::select(-id) %>%
        dplyr::rename_with(.fn =  ~gsub("pr", "prcp_", .x),
                           .cols = dplyr::starts_with("pr")) %>%
        dplyr::rename_with(.fn =  ~gsub("tmmn", "tmmn_", .x),
                           .cols = dplyr::starts_with("tmmn")) %>%
        dplyr::rename_with(.fn =  ~gsub("tmmx", "tmmx_", .x),
                           .cols = dplyr::starts_with("tmmx")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("tmmn"),
                                    .fns = ~.x/10),
                      dplyr::across(.cols = dplyr::starts_with("tmmx"),
                                    .fns = ~.x/10)) %>%
        sf::st_drop_geometry()


    # Annotate with yearly surface water occurrence --------------------------------

    # Number of pixels with water each year
    band <- "waterClass"
    stackCollection <- ABAP::EEcollectionToMultiband(collection = "JRC/GSW1_3/YearlyHistory",
                                                     dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                     band = band,
                                                     group_type = "year",
                                                     groups = config$years,
                                                     unmask = FALSE)

    out <- ABAP::addVarEEimage(ee_pentads, stackCollection, "count")

    # Fix covariates
    out <- out %>%
        dplyr::select(Name, dplyr::everything()) %>%
        dplyr::select(-id) %>%
        dplyr::rename_with(~gsub("X", "watext_", .x), .cols = dplyr::starts_with("X")) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")

    # Recurrence of pixels with water each year
    out <- ABAP::addVarEEimage(ee_pentads, stackCollection, "mean")

    # Fix covariates
    out <- out %>%
        dplyr::select(Name, dplyr::everything()) %>%
        dplyr::select(-id) %>%
        dplyr::rename_with(~gsub("X", "watrec_", .x), .cols = dplyr::starts_with("X")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("wat"),
                                    .fns = ~tidyr::replace_na(.x, 0))) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")


    # Annotate with overall surface water occurrence --------------------------

    out <- ABAP::addVarEEimage(ee_pentads = ee_pentads,
                               image = "JRC/GSW1_3/GlobalSurfaceWater",
                               reducer = "mean",
                               bands = "occurrence",
                               unmask = TRUE)

    # Fix covariates
    out <- out %>%
        dplyr::select(Name, dplyr::everything()) %>%
        dplyr::select(-id) %>%
        dplyr::rename(water_occur = occurrence) %>%
        dplyr::rename_with(.fn =  ~gsub("water_occur", "watocc_ever", .x),
                           .cols = dplyr::starts_with("water_occur")) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")


    # Annotate with yearly NDVI -----------------------------------------------

    # Find mean water presence for each pixel and year
    band <- "NDVI"
    stackCollection <- ABAP::EEcollectionToMultiband(collection = "MODIS/006/MOD13A2",
                                                     dates = paste0(config$year_range + c(0,1), "-01-01"),
                                                     band = band,
                                                     group_type = "year",
                                                     groups = config$years,
                                                     reducer = "mean",
                                                     unmask = FALSE)

    # Find mean (mean) water presence for each pentad and year
    out <- ABAP::addVarEEimage(ee_pentads, stackCollection, "mean")

    # Fix covariates
    out <- out %>%
        dplyr::select(Name, dplyr::everything()) %>%
        dplyr::select(-id) %>%
        dplyr::rename_with(~gsub("X", "ndvi_", .x), .cols = dplyr::starts_with("X")) %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("ndvi"),
                                    .fns = ~.x/10000)) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        dplyr::left_join(out, by = "Name")


    # Add geometries -----------------------------------------------------------

    sitedata <- sitedata %>%
        dplyr::left_join(dplyr::select(pentads_sa, Name), by = "Name")

    # Add centroids
    # sf::st_agr(sitedata) = "constant"
    cc <- sf::st_centroid(sitedata) %>%
        sf::st_coordinates() %>%
        as.data.frame()

    sitedata <- sitedata %>%
        dplyr::mutate(lon = cc$X,
                      lat = cc$Y)


    # Annotate with distance to coast -----------------------------------------

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

    out <- ABAP::addVarEEimage(ee_pentads = ee_pentads,
                               image = "MERIT/DEM/v1_0_3",
                               reducer = "mean",
                               bands = "dem",
                               unmask = FALSE)

    # Fix covariates
    out <- out %>%
        dplyr::select(Name, elev = dem) %>%
        sf::st_drop_geometry()

    sitedata <- sitedata %>%
        sf::st_drop_geometry() %>%
        dplyr::left_join(out, by = "Name")

    utils::write.csv(sitedata,
              file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv")),
              row.names = FALSE)

    return(sitedata)

}
