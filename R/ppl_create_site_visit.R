#' Create site and visit occupancy files
#'
#' @description This function prepares site and visit occupancy data to fit an
#' occupancy model from ABAP data. It has two parts: the first part downloads ABAP
#' data and annotates them with covariates from Google Earth Engine using the
#' functions \code{\link{prepGEESiteData()}} and \code{\link{prepGEEVisitData()}},
#' the second part uses the function \code{\link{createOccuData}} to format the
#' data. Other changes might be necessary to fit into different occupancy modelling
#' packages
#' @inheritParams ppl_run_pipe_dst1
#' @param force_gee_dwld Whether covariates from Google Earth Engine should be
#' downloaded, even if a file with covariates is already present on disk.
#' Defaults to FALSE.
#' @param save_occu_data If TRUE data is saved to disc, but see 'overwrite_occu_data'.
#' @param overwrite_occu_data A character vector with the data that should be overwritten in
#' case it is already present on disc. It can be any combination of
#' c("site", "visit", "det"), with site referring to site data, visit to visit data
#' and det to detection data. Site and visit data are typically common for
#' multiple species and we might not want to save it all the time.
#'
#' @return The first part of the function creates two data frames (in .csv format)
#' that will be saved to disk: GEE annotated ABAP site data and GEE annotated ABAP
#' visit data. The second part of the functions creates three data frames that will
#' be saved to disk: site, visit and species detection data frames,
#' all in .csv format.
#' @export
#'
#' @examples
ppl_create_site_visit <- function(sp_code, force_gee_dwld = FALSE,
                                  save_occu_data = TRUE,
                                  overwrite_occu_data = NULL, config, ...){

    varargs <- list(...)


    # Download ABAP data and annotate with Google Earth Engine ----------------

    sitefile <- file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv"))
    visitfile <- file.path(config$out_dir, paste0("visit_dat_sa_gee_", config$years_ch, ".csv"))

    # Download from GEE if file doesn't exit
    if("force_gee_dwld" %in% names(varargs)){
        force_gee_dwld <- varargs$force_gee_dwld
    }

    if(!file.exists(sitefile) | force_gee_dwld){
        prepGEESiteData(config, monitor = varargs$monitor_gee)
    }

    if(!file.exists(visitfile) | force_gee_dwld){
        prepGEEVisitData(config, monitor = varargs$monitor_gee)
    }

    # Load data
    sitedata <- utils::read.csv(sitefile)
    visitdata <- utils::read.csv(visitfile)


    # Format to occu ---------------------------------------------------------

    occudata <- BIRDIE::createOccuData(sp_code = sp_code,
                                       years = config$years,
                                       site_data = sitedata,
                                       visit_data = visitdata,
                                       config = config,
                                       force_abap_dwld = varargs$force_abap_dwld)


    # Subset sites that have been visited in the period
    rm(sitedata, visitdata)
    site_data <- occudata$site
    visit_data <- occudata$visit
    rm(occudata)
    gc()

    # Keep only pentads that appear in visit data
    site_data <- site_data %>%
        dplyr::filter(Pentad %in% unique(visit_data$Pentad))

    # Check that all pentads-years in visit data are also in site data
    miss_pentads <- visit_data %>%
        dplyr::as_tibble() %>%
        dplyr::distinct(Pentad, year) %>%
        dplyr::anti_join(
            site_data %>%
                dplyr::as_tibble() %>%
                dplyr::distinct(Pentad, year),
            by = c("Pentad", "year"))

    if(nrow(miss_pentads) > 0){
        log_name <- file.path(config$out_dir, "reports", paste0("pentads_in_visit_not_site_", config$years_ch, ".csv"))
        utils::write.csv(miss_pentads, log_name, row.names = FALSE)
    }


    # Remove visit data from missing pentads
    visit_data <- visit_data %>%
        dplyr::anti_join(miss_pentads, by = c("Pentad", "year"))

    if(save_occu_data){

        visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
        sitefile <- file.path(config$out_dir, paste0("occu_site_dat_sa_", config$years_ch, ".csv"))
        detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

        if((!file.exists(visitfile)) | (file.exists(visitfile) & ("visit" %in% overwrite_occu_data))){
            visit_data %>%
                dplyr::select(-obs) %>%
                utils::write.csv(visitfile, row.names = FALSE)
        } else {
            warning("Visit file not saved because 'visit' not in overwrite_occu_data")
        }

        if((!file.exists(sitefile)) | (file.exists(sitefile) & ("site" %in% overwrite_occu_data))){
            site_data %>%
                utils::write.csv(sitefile, row.names = FALSE)
        } else {
            warning("Site file not saved because 'site' not in overwrite_occu_data")
        }

        if((!file.exists(detfile)) | (file.exists(detfile) & ("det" %in% overwrite_occu_data))){
            visit_data %>%
                dplyr::select(CardNo, StartDate, Pentad, year, obs) %>%
                utils::write.csv(detfile, row.names = FALSE)
        } else {
            warning("Detection file not saved because 'det' not in overwrite_occu_data")
        }

    }

}
