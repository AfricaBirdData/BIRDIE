#' Create site and visit occupancy files
#'
#' @inheritParams ppl_run_pipe_distr
#'
#' @return
#' @export
#'
#' @examples
ppl_create_site_visit <- function(sp_code, year, config, ...){

    varargs <- list(...)

    # Data file name
    if(year < 2020){
        datafile <- file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")
    } else {
        yy <- substring(as.character(config$year), 3, 4)
        ff <- paste0("site_dat_sa_gee_", yy, ".rds")
        datafile <- file.path(config$data_dir, ff)
    }

    # Load data and subset years
    sitedata <- readRDS(datafile) %>%
        sf::st_drop_geometry() %>%
        dplyr::select(Pentad = Name, lon, lat, watocc_ever, dist_coast, dplyr::ends_with(match = as.character(config$years))) %>%
        tidyr::drop_na()   # I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE

    visitdata <- readRDS(file.path(config$data_dir, "visit_dat_sa_gee_08_19.rds")) %>%
        dplyr::filter(year %in% config$years)


    # Format to occuR ---------------------------------------------------------

    occuRdata <- BIRDIE::prepDataOccuR(spp_code = sp_code,
                                       years = config$years,
                                       site_data = sitedata,
                                       visit_data = visitdata,
                                       config = config,
                                       download = varargs$download_from_abap,
                                       save = varargs$save_occu_data,
                                       overwrite = varargs$overwrite_occu_data)

    return(occuRdata)

}
