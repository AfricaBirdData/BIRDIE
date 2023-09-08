#' Create site and visit occupancy files
#'
#' @description This function prepares site and visit occupancy data to fit an
#' occupancy model from ABAP data. It has two parts: the first part downloads ABAP
#' data and annotates them with covariates from Google Earth Engine using the
#' functions \code{\link{prepGEESiteData}} and \code{\link{prepGEEVisitData}},
#' the second part uses the function \code{\link{createOccuData}} to format the
#' data.
#' @inheritParams ppl_run_pipe_dst1
#' @param force_gee_dwld Whether covariates from Google Earth Engine should be
#' downloaded, even if a file with covariates is already present on disk.
#' Defaults to FALSE.
#' @param force_site_visit Whether site and visit data should be prepared even if
#' visit and site data files are already on disk. Defaults to FALSE
#' @param force_abap_dwld Indicates whether ABAP data must be downloaded for the species
#' and years indicated by 'sp_code' and 'years'. If TRUE, data will be downloaded
#' from ABAP once per session and cached in a temp file. After this the cached
#' file will be used, unless download is set to FALSE, in which case data will
#' be downloaded regardless of the cached file.
#' @param monitor_gee if TRUE (default) periodic messages of the state of the downloads
#' from GEE will be printed on screen.
#'
#' @return The first part of the function creates two data frames (in .csv format)
#' that will be saved to disk: GEE annotated ABAP site data and GEE annotated ABAP
#' visit data. The second part of the functions creates three data frames that will
#' be saved to disk: site, visit and species detection data frames,
#' all in .csv format.
#' @export
#'
#' @examples
ppl_create_site_visit <- function(config,
                                  sp_code,
                                  force_gee_dwld = FALSE,
                                  force_site_visit = FALSE,
                                  force_abap_dwld = FALSE,
                                  monitor_gee = TRUE){

    # varargs <- list(...)

    # Should site and visit data preparation be forced?
    visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
    sitefile <- file.path(config$out_dir, paste0("occu_site_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

    if(force_site_visit){
        prep_data <- TRUE
    } else if(!file.exists(sitefile) | !file.exists(visitfile)){
        prep_data <- TRUE
    } else {
        prep_data <- FALSE
    }


    if(prep_data){

        # Download ABAP data and annotate with Google Earth Engine ----------------

        geesitefile <- file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv"))
        geevisitfile <- file.path(config$out_dir, paste0("visit_dat_sa_gee_", config$years_ch, ".csv"))

        if(!file.exists(geesitefile) | !file.exists(geevisitfile) | force_gee_dwld){
            # Initialize Earth Engine
            rgee::ee_check()
            rgee::ee_Initialize(drive = TRUE)
        }

        # Download from GEE if file doesn't exit
        if(!file.exists(geesitefile) | force_gee_dwld){

            # The following code is HARD CODED. In the future we might want to run
            # the pipeline for areas other than SA and we might need to upload
            # new pentads to GEE (ie. set upload_asset = TRUE)
            pentads_sa <- ABAP::getRegionPentads(.region_type = "country",
                                                 .region = "South Africa")

            geesitedata <- prepGEESiteData(config = config,
                                           pentads = pentads_sa,
                                           asset_id = "birdie_pentads_sa",  # HARD CODED
                                           upload_asset = TRUE,
                                           monitor = monitor_gee)

            utils::write.csv(geesitedata, geesitefile, row.names = FALSE)

            message(paste("Site data with GEE covts saved at", geesitefile))

        } else {
            # Load data
            geesitedata <- utils::read.csv(geesitefile)
        }

        if(!file.exists(geevisitfile) | force_gee_dwld){

            # Prepare visit data to annotate with GEE environmental data
            pentads_sa <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa") # HARD CODED

            years <- config$years

            visits <- data.frame()

            for(i in seq_along(years)){

                # Download SABAP data for any species

                year <- years[i]

                visit <- ABAP::getAbapData(.spp_code = 6,  # this species should not matter - visits are the same for all species
                                           .region_type = "country",
                                           .region = "South Africa",
                                           .years = year)

                visits <- dplyr::bind_rows(visits, visit)

            }

            # Make spatial object and select relevant columns
            visits <- visits %>%
                dplyr::left_join(pentads_sa,
                                 by = c("Pentad" = "Name")) %>%
                sf::st_sf() %>%
                dplyr::filter(!sf::st_is_empty(.)) %>%     # Remove rows without geometry
                dplyr::mutate(year = lubridate::year(StartDate),
                              Date = as.character(StartDate)) %>%   # GEE doesn't like dates
                dplyr::select(year, CardNo, StartDate, Date, Pentad, TotalHours, ObserverNo)

            geevisitdata <- prepGEEVisitData(config = config,
                                             visits = visits,
                                             asset_id = "birdie_visits",
                                             upload_asset = TRUE,
                                             monitor = monitor_gee)

            utils::write.csv(geevisitdata, geevisitfile, row.names = FALSE)

            message(paste("Visits data with GEE covts saved at", geevisitfile))

        } else {
            # Load data
            geevisitdata <- utils::read.csv(geevisitfile)
        }


        # Prepare site and visit data ---------------------------------------------

        occudata <- BIRDIE::createOccuData(config = config,
                                           sp_code = sp_code,
                                           years = config$years,
                                           site_data = geesitedata,
                                           visit_data = geevisitdata)

        # Clean environment
        rm(geesitedata, geevisitdata)
        site_data <- occudata$site
        visit_data <- occudata$visit
        rm(occudata)
        gc()

    } else {

        site_data <- read.csv(sitefile, check.names = FALSE)
        visit_data <- read.csv(visitfile, check.names = FALSE)

    }

    # Download detection data -------------------------------------------------

    # Cache file name
    cachefile <- file.path(tempdir(), paste(c(sp_code, config$years, ".rds"), collapse = "_"))

    if(!file.exists(cachefile) | force_abap_dwld){

        # Download species detection
        message("Downloading from SABAP")

        sp_detect <- ABAP::getAbapData(.spp_code = sp_code,
                                       .region_type = "country",
                                       .region = "South Africa",
                                       .years = config$years)

        # Save to cache
        saveRDS(sp_detect, cachefile)

    } else {

        message("Using cached file")
        sp_detect <- readRDS(cachefile)

    }

    # Add detections to visit data
    visit_data <- visit_data %>%
        dplyr::mutate(StartDate = as.Date(StartDate)) %>%
        dplyr::left_join(sp_detect %>%
                             dplyr::select(CardNo, StartDate, Pentad, obs = Spp) %>%
                             dplyr::mutate(obs = ifelse(obs == "-", 0, 1)),
                         by = c("CardNo", "StartDate", "Pentad"))

    # Remove NA values in detection data
    if(any(is.na(visit_data$obs))){
        warning("NA found in detection data")
        visit_data <- visit_data %>%
            dplyr::filter(!is.na(obs))
    }

    return(list(site = site_data,
                visit = visit_data))

}
