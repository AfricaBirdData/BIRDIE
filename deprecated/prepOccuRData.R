#' Prepare occuR data for model fitting
#'
#' @inheritParams ppl_run_pipe_distr
#'
#' @return
#' @export
#'
#' @examples
prepOccuRData <- function(sp_code, year, config, ...){

    varargs <- list(...)

    # Prepare data ------------------------------------------------------------

    # File names
    visitfile <- file.path(config$out_dir, paste0("occur_visit_dat_sa_", config$years_ch, ".csv"))
    sitefile <- file.path(config$out_dir, paste0("occur_site_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occur_det_dat_sa_", config$years_ch, ".csv"))

    # Read in site and visit data
    site_data <- read.csv(sitefile)
    visit_data <- read.csv(visitfile)
    det_data <- read.csv(detfile)

    # Stop if there are no detections
    if(!1 %in% unique(det_data$obs)){
        warning(paste("No detection of species", sp_code))
        return(1)
    }

    # Or species detected in too few Pentads
    n_pentads <- det_data %>%
        dplyr::count(Pentad, obs) %>%
        dplyr::filter(obs == 1) %>%
        nrow()

    if(n_pentads < 5){
        warning(paste("Species", sp_code, "detected in less than 5 pentads"))
        return(2)
    }

    # Add detection data to visits
    visit_data <- visit_data %>%
        dplyr::left_join(det_data, by = c("Pentad", "year", "visit"))

    site_data <- data.table::as.data.table(site_data)
    visit_data <- data.table::as.data.table(visit_data)


    # Scaling -----------------------------------------------------------------

    occuRdata <- scaleOccuRVars(site_data, visit_data,
                                scale_vars = varargs$scale_vars_occur)

    return(occuRdata)

}
