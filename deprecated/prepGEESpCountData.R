#' Prepare Google Earth Engine data for CWAC species counts
#'
#' @param counts A dataframe with counts to annotate with covariates. This data
#' frame must contain at least the lat-lon coordinates of the CWAC site the
#' counts were obtained at. The variables must be named "X" and "Y".
#' @param sp_code SAFRING reference number of the species we want to analyze.
#' @param catchment An sf object with the polygons defining the catchments
#' of the CWAC sites in `counts`. Which site belongs to which catchment will be
#' determined by a join in GEE.
#' @param config A list with pipeline configuration parameters.
#' See \link{configPipeline}
#' @param monitor Logical. If TRUE (default) monitoring printed messages produced
#' by `rgee` will displayed. If FALSE, only high-level messages will be displayed.
#' @param force_gee If TRUE annotation with GEE info will occur even if a
#' dataset with covariates is already present on disk.
#'
#' @return
#' @export
#'
#' @examples
prepGEESpCountData <- function(counts, sp_code, catchment, config, monitor = TRUE,
                               force_gee = FALSE){

    rgee::ee_check()
    rgee::ee_Initialize()

    # Set species code and output file name
    outfile <- setSpOutFilePath("abu_gee_data", config, sp_code, ".csv")

    if(!force_gee & file.exists(outfile)){
        stop("File with covariates already on disk. Set force_gee = TRUE to overwrite")
    }

    # Set a name for the asset
    eeCatchm_id <- file.path(rgee::ee_get_assethome(), 'quin_catchm')

    # Find nearest catchment to site counts. Nearest because some sites are on the
    # coast potentially falling offshore
    counts <- counts %>%
        dplyr::mutate(id_count = dplyr::row_number()) %>%
        sf::st_as_sf(coords = c("X", "Y"), dim = "XY", crs = sf::st_crs(4326))

    int_index <- counts %>%
        sf::st_nearest_feature(catchment)

    counts <- counts %>%
        dplyr::mutate(UNIT_ID = catchment$UNIT_ID[int_index])


    # Upload to GEE -----------------------------------------------------------

    # Set a name for the asset
    eeCounts_id <- file.path(rgee::ee_get_assethome(), 'birdie_cwac_counts')

    # Upload to EE (if not done already). Delete first
    assets <- rgee::ee_manage_assetlist(rgee::ee_get_assethome())

    if(eeCounts_id %in% assets$ID){
        rgee::ee_manage_delete(eeCounts_id, quiet = FALSE, strict = TRUE)
    }

    counts %>%
        dplyr::rename(Date = StartDate) %>%
        dplyr::mutate(Date = lubridate::floor_date(Date, "month")) %>%
        dplyr::mutate(Date = as.character(Date)) %>%
        dplyr::select(id_count, Date, UNIT_ID) %>%
        ABDtools::uploadFeaturesToEE(asset_id = eeCounts_id,
                                     load = FALSE,
                                     monitor = monitor)

    eeCounts <- rgee::ee$FeatureCollection(eeCounts_id)

    eeCatchment <- rgee::ee$FeatureCollection(eeCatchm_id)

    # Use an equals filter to specify how the collections match.
    eeFilter = rgee::ee$Filter$equals(
        leftField = 'UNIT_ID',
        rightField = 'UNIT_ID')

    # Define the join.
    saveAllJoin <- rgee::ee$Join$saveAll(
        matchesKey = 'UNIT_ID',
        ordering = 'system:time_start',
        ascending = TRUE)

    # Apply the join
    joinResults <- saveAllJoin$apply(eeCounts, eeCatchment, eeFilter)

    # Extract polygons from join results
    new_pols <- joinResults$map(function(a){

        pp <- rgee::ee$Feature(rgee::ee$List(a$get("UNIT_ID"))$get(0))
        ll <- rgee::ee$List(a)

        return(
            pp$copyProperties(ll, list("id_count", "Date"))
        )

    })


    # Annotate with TerraClimate ----------------------------------------------

    message("Annotating CWAC with TerraClimate")

    # Define bands
    bands <- c("pr", "tmmn", "tmmx", "pdsi")

    # Define function
    f <- function(band){

        # Annotate with GEE TerraClimate
        visit_env <- ABDtools::addVarEEclosestImage(ee_feats = new_pols,
                                                    collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                    reducer = "mean",
                                                    maxdiff = 15,
                                                    bands = band,
                                                    monitor = monitor)

        # Fix names and variables
        visit_env <- visit_env %>%
            dplyr::select(id_count, dplyr::all_of(paste0(band, "_mean"))) %>%
            sf::st_drop_geometry()

        return(visit_env)

    }

    # Annotate
    visit_vars <- purrr::map(bands, ~f(.x))

    visit_vars <- visit_vars %>%
        dplyr::bind_cols(.name_repair = "universal") %>%
        dplyr::rename(id_count = id_count...1) %>%
        dplyr::select(id_count, dplyr::all_of(paste0(bands, "_mean")))

    # bind
    counts_vars <- counts %>%
        dplyr::left_join(visit_vars,
                         by = "id_count")

    # Fix covariates
    counts_vars <- visit_vars %>%
        dplyr::select(id_count, dplyr::everything()) %>%
        dplyr::rename_with(.fn =  ~gsub("pr", "prcp", .x),
                           .cols = dplyr::starts_with("pr")) %>%
        dplyr::rename_with(.fn =  ~gsub("tmmn", "tmmn", .x),
                           .cols = dplyr::starts_with("tmmn")) %>%
        dplyr::rename_with(.fn =  ~gsub("tmmx", "tmmx", .x),
                           .cols = dplyr::starts_with("tmmx")) %>%
        dplyr::mutate(across(.cols = dplyr::starts_with("tmmn"),
                             .fns = ~.x/10),
                      across(.cols = dplyr::starts_with("tmmx"),
                             .fns = ~.x/10))


    # Annotate with water occurrence ------------------------------------------

    message("Annotating CWAC with surface water")

    visit_water <- vector("list", length = 2)

    # Number of pixels with water each year
    visit_water[[1]] <- ABDtools::addVarEEclosestImage(ee_feats = new_pols,
                                                       collection = "JRC/GSW1_3/YearlyHistory",
                                                       reducer = "count",
                                                       maxdiff = 1000,
                                                       bands = "waterClass",
                                                       monitor = monitor)

    # Fix names and variables
    visit_water[[1]] <- visit_water[[1]] %>%
        dplyr::rename_with(~gsub("waterClass", "watext", .x), .cols = dplyr::starts_with("waterClass")) %>%
        dplyr::select(id_count, dplyr::all_of("watext_count")) %>%
        sf::st_drop_geometry()

    # Recurrence of pixels with water each year
    visit_water[[2]] <- ABDtools::addVarEEclosestImage(ee_feats = new_pols,
                                                       collection = "JRC/GSW1_3/YearlyHistory",
                                                       reducer = "mean",
                                                       maxdiff = 1000,
                                                       bands = "waterClass",
                                                       monitor = monitor)

    # Fix names and variables
    visit_water[[2]] <- visit_water[[2]] %>%
        dplyr::rename_with(~gsub("waterClass", "watrec", .x), .cols = dplyr::starts_with("waterClass")) %>%
        dplyr::select(id_count, dplyr::all_of("watrec_mean")) %>%
        sf::st_drop_geometry()

    # bind
    visit_vars <- visit_water[[1]] %>%
        dplyr::left_join(visit_water[[2]],
                         by = "id_count")

    # join with previous variables and save
    counts_vars <- counts_vars %>%
        dplyr::left_join(visit_vars,
                         by = "id_count")


    # Return data with covariates -----------------------------------------------

    # Correct
    counts <- counts %>%
        sf::st_drop_geometry() %>%
        dplyr::left_join(counts_vars, by = "id_count")

    return(counts)

}
