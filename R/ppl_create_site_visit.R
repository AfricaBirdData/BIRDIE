#' Create site and visit occupancy files
#'
#' @description This function prepares site and visit occupancy data to fit an
#' occupancy model from ABAP data. It has two parts: the first part downloads ABAP
#' data and annotates them with covariates from Google Earth Engine using the
#' functions \code{\link{prepGEESiteData()}} and \code{\link{prepGEEVisitData()}},
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
#'
#' @return The first part of the function creates two data frames (in .csv format)
#' that will be saved to disk: GEE annotated ABAP site data and GEE annotated ABAP
#' visit data. The second part of the functions creates three data frames that will
#' be saved to disk: site, visit and species detection data frames,
#' all in .csv format.
#' @export
#'
#' @examples
ppl_create_site_visit <- function(config, sp_code,
                                  force_gee_dwld = FALSE,
                                  force_site_visit = FALSE,
                                  force_abap_dwld = FALSE){

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

        # Download from GEE if file doesn't exit
        if(!file.exists(geesitefile) | force_gee_dwld){
            prepGEESiteData(config, monitor = varargs$monitor_gee)
        }

        if(!file.exists(geevisitfile) | force_gee_dwld){
            prepGEEVisitData(config, monitor = varargs$monitor_gee)
        }

        # Load data
        geesitedata <- utils::read.csv(geesitefile)
        geevisitdata <- utils::read.csv(geevisitfile)


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

        site_data <- read.csv(sitefile)
        visit_data <- read.csv(visitfile)

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
