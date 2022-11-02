#' Create site and visit occupancy files
#'
#' @description This function prepares site and visit occupancy data to fit an
#' occuR model from ABAP data. It has two parts: the first part downloads ABAP
#' data and annotates them with covariates from Google Earth Engine using the
#' functions \code{\link{prepGEESiteData()}} and \code{\link{prepGEEVisitData()}},
#' the second part uses the function \code{\link{createOccuRData}} to format the
#' data for the \code{occuR} package.
#' @inheritParams ppl_run_pipe_dst1
#' @param force_gee_dwld Whether covariates from Google Earth Engine should be
#' downloaded, even if a file with covariates is already present on disk.
#' Defaults to FALSE.
#'
#' @return The first part of the function creates two data frames (in .csv format)
#' that will be saved to disk: GEE annotated ABAP site data and GEE annotated ABAP
#' visit data. The second part of the functions creates three data frames that will
#' be saved to disk: occuR-formatted site, visit and species detection data frames,
#' all in .csv format.
#' @export
#'
#' @examples
ppl_create_site_visit <- function(sp_code, year, force_gee_dwld = FALSE,
                                  config, ...){

    varargs <- list(...)

    # Site and visit data file names
    if(year < 2020){
        sitefile <- file.path(config$out_dir, "site_dat_sa_gee_08_19.csv")
        visitfile <- file.path(config$out_dir, "visit_dat_sa_gee_08_19.csv")
    } else {
        sitefile <- file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv"))
        visitfile <- file.path(config$out_dir, paste0("visit_dat_sa_gee_", config$years_ch, ".csv"))
    }

    # Download from GEE if file doesn't exit
    if("force_gee_dwld" %in% names(varargs)){
        force_gee_dwld <- varargs$force_gee_dwld
    }

    if(!file.exists(sitefile) | force_gee_dwld){
        prepGEESiteData(config, monitor = varargs$monitor)
    }

    if(!file.exists(visitfile) | force_gee_dwld){
        prepGEEVisitData(config, monitor = varargs$monitor)
    }

    # Load data and subset years
    sitedata <- utils::read.csv(sitefile) %>%
        dplyr::select(Pentad = Name, lon, lat, watocc_ever, dist_coast, dplyr::ends_with(match = as.character(config$years))) %>%
        tidyr::drop_na()   # I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE

    visitdata <- utils::read.csv(visitfile) %>%
        dplyr::filter(year %in% config$years)


    # Format to occuR ---------------------------------------------------------

    occuRdata <- BIRDIE::createOccuRData(sp_code = sp_code,
                                         years = config$years,
                                         site_data = sitedata,
                                         visit_data = visitdata,
                                         config = config,
                                         force_abap_dwld = varargs$force_abap_dwld,
                                         save_occu_data = varargs$save_occu_data,
                                         overwrite_occu_data = varargs$overwrite_occu_data)

    return(occuRdata)

}
