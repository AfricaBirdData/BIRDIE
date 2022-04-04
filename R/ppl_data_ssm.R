#' Prepare covariate data for abundance pipeline
#'
#' @inheritParams ppl_run_pipe_abu1
#'
#' @return
#' @export
#'
#' @examples
ppl_data_ssm <- function(site, year, config, ...){

    rgee::ee_check()
    rgee::ee_Initialize(drive = TRUE)


    # Prepare site data -------------------------------------------------------

    # List all sites in South Africa
    sites <- CWAC::listCwacSites(.region_type = "country", .region = "South Africa")

    # Extract site of interest
    site_spt <- sites %>%
        dplyr::filter(LocationCode == site) %>%
        dplyr::select(LocationCode, LocationName, X, Y) %>%
        sf::st_as_sf(coords = c("X", "Y"), dim = "XY", crs = sf::st_crs(4326))

    # Read in catchment data
    catchmt <- sf::read_sf(file.path(config$data_dir, "catchmt_4.shp")) # THIS IS HARD CODED

    # Extract catchment
    site_ctm <- catchmt[unlist(sf::st_intersects(site_spt, catchmt)),]


    # Get counts for the site -------------------------------------------------

    # Download counts for the site
    counts <- CWAC::getCwacSiteCounts(site)

    # Add Doug's data if site is Du Toit's Pan
    if(site == 28462448){
        counts <- rbind(counts,
                        utils::read.csv(file.path(config$data_dir, "28462448_data_2022_doug.csv"))) %>%
            dplyr::arrange(StartDate)
    }

    # Set analysis period as the period from first count to two years after last count
    year_rg <- range(counts$Year)

    yy_to_present <- lubridate::year(Sys.time()) - year_rg[2] - 1

    if(yy_to_present > 2){
        year_rg[2] <- year_rg[2] + 2
    }

    counts <- counts %>%
        dplyr::filter(Year >= year_rg[1], Year <= year_rg[2])

    # Add missing surveys
    counts <- CWAC::addMissingCwacCounts(counts, years = year_rg[1]:year_rg[2])

    # Give missing surveys a date based on the dates from other surveys
    month_summer <- counts %>%
        dplyr::mutate(month = lubridate::month(StartDate)) %>%
        dplyr::filter(Season == "S", !is.na(month)) %>%
        dplyr::count(month) %>%
        dplyr::filter(n == max(n)) %>%
        dplyr::pull(month)

    month_winter <- counts %>%
        dplyr::mutate(month = lubridate::month(StartDate)) %>%
        dplyr::filter(Season == "W", !is.na(month)) %>%
        dplyr::count(month) %>%
        dplyr::filter(n == max(n)) %>%
        dplyr::pull(month)

    counts <- counts %>%
        dplyr::mutate(StartDate = dplyr::case_when(is.na(StartDate) & Season == "S" ~ as.Date(paste(Year, month_summer, "01", sep = "-")),
                                                   is.na(StartDate) & Season == "W" ~ as.Date(paste(Year, month_winter, "01", sep = "-")),
                                                   TRUE ~ StartDate) ) %>%
        dplyr::arrange(StartDate, Season)

    # Extract visit dates, add geometry and id
    visits <- counts %>%
        dplyr::distinct(Card, StartDate) %>%
        dplyr::mutate(id = dplyr::row_number(),
                      geometry = site_ctm$geometry) %>%
        sf::st_sf()

    # Create output name
    outfile <- file.path(config$data_dir, paste(site, year, "visit_covts.rds", sep = "_"))


    # Upload to GEE -----------------------------------------------------------

    # Set a name for the asset
    eeid <- sprintf("%s/%s", rgee::ee_get_assethome(), 'cwac_visits')

    # Upload to EE (if not done already)
    visits %>%
        dplyr::rename(Date = StartDate) %>%
        dplyr::mutate(Date = lubridate::floor_date(Date, "month")) %>%
        dplyr::mutate(Date = as.character(Date)) %>%
        ABAP::uploadPentadsToEE(asset_id = eeid,
                                load = FALSE)


    # Annotate with TerraClimate ----------------------------------------------

    # Load the remote data asset
    ee_visit <- rgee::ee$FeatureCollection(eeid)

    # Define bands
    bands <- c("pr", "tmmn", "tmmx")

    # Define function
    f <- function(band){

        # Annotate with GEE TerraClimate
        visit_env <- ABAP::addVarEEclosestImage(ee_pentads = ee_visit,
                                                collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                                reducer = "mean",
                                                maxdiff = 15,
                                                bands = band)

        # Fix names and variables
        visit_env <- visit_env %>%
            dplyr::rename_with(~gsub("val", band, .x), .cols = dplyr::starts_with("val")) %>%
            dplyr::select(id, dplyr::all_of(band)) %>%
            sf::st_drop_geometry()

        return(visit_env)

    }

    # Annotate
    visit_vars <- purrr::map(bands, ~f(.x))

    # bind
    visit_vars <- visits %>%
        sf::st_drop_geometry() %>%
        dplyr::left_join(dplyr::bind_cols(visit_vars) %>%
                             dplyr::rename(id = id...1) %>%
                             dplyr::select(id, dplyr::all_of(bands)),
                         by = "id")

    # Fix covariates
    visit_vars <- visit_vars %>%
        dplyr::select(id, dplyr::everything()) %>%
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

    saveRDS(visit_vars, outfile)


    # Annotate with water occurrence ------------------------------------------

    # Load the remote data asset
    ee_visit <- rgee::ee$FeatureCollection(eeid)

    visit_water <- vector("list", length = 2)

    # Number of pixels with water each year
    visit_water[[1]] <- ABAP::addVarEEclosestImage(ee_pentads = ee_visit,
                                                   collection = "JRC/GSW1_3/YearlyHistory",
                                                   reducer = "count",
                                                   maxdiff = 1000,
                                                   bands = "waterClass")

    # Fix names and variables
    visit_water[[1]] <- visit_water[[1]] %>%
        dplyr::rename_with(~gsub("val", "watext", .x), .cols = dplyr::starts_with("val")) %>%
        dplyr::select(id, dplyr::all_of("watext")) %>%
        sf::st_drop_geometry()

    # Recurrence of pixels with water each year
    visit_water[[2]] <- ABAP::addVarEEclosestImage(ee_pentads = ee_visit,
                                                   collection = "JRC/GSW1_3/YearlyHistory",
                                                   reducer = "mean",
                                                   maxdiff = 1000,
                                                   bands = "waterClass")

    # Fix names and variables
    visit_water[[2]] <- visit_water[[2]] %>%
        dplyr::rename_with(~gsub("val", "watrec", .x), .cols = dplyr::starts_with("val")) %>%
        dplyr::select(id, dplyr::all_of("watrec")) %>%
        sf::st_drop_geometry()

    # bind
    visit_vars <- visit_water[[1]] %>%
        dplyr::left_join(dplyr::bind_cols(visit_water[[2]]),
                         by = "id")

    # join with previous variables and save
    readRDS(outfile) %>%
        dplyr::left_join(visit_vars, by = "id") %>%
        saveRDS(outfile)


    # Save data with covariates -----------------------------------------------

    visit_vars <- readRDS(outfile)

    counts %>%
        dplyr::left_join(dplyr::select(visit_vars, -id), by = c("Card", "StartDate")) %>%
        saveRDS(outfile)

}
