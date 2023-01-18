#' Fit occupancy model
#'
#' @inheritParams ppl_run_pipe_dst1
#' @param year_sel Year in data to run model for.
#' @param spatial Whether a spatial model should be fit. Defaults to FALSE.
#'
#' @return It can return either a number correspoding to the status of the fitting
#' process: 1 - No detections, 2 - too few detections, 3 - model fitting failed
#' (inherited from \code{\link{fitSpOccu}})
#' @export
#'
#' @examples
ppl_fit_occu_model <- function(sp_code, year_sel, config, spatial = FALSE, ...){

    varargs <- list(...)

    # Prepare data ------------------------------------------------------------

    # File names
    visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
    sitefile <- file.path(config$out_dir, paste0("occu_site_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

    # Read in site and visit data
    site_data <- utils::read.csv(sitefile, check.names = FALSE)
    visit_data <- utils::read.csv(visitfile, check.names = FALSE)
    det_data <- utils::read.csv(detfile, check.names = FALSE)

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

    # Add detection info to visit data (we need Spp variable for functions in the ABAP package)
    visit_data <- visit_data %>%
        dplyr::left_join(det_data,
                         by = c("CardNo", "StartDate", "year", "Pentad")) %>%
        dplyr::mutate(Spp = ifelse(obs == 0, "-", sp_code))

    # Subset data sets
    site_data_year <- site_data %>%
        dplyr::filter(year == year_sel)

    visit_data_year <- visit_data %>%
        dplyr::filter(year == year_sel)

    # Get spatial sites if necessary
    if(spatial){
        sp_sites <- ABAP::getRegionPentads("country", "South Africa") %>%
            dplyr::filter(Name %in% unique(site_data_year$Name)) %>%
            dplyr::select(Name, pentad)
    } else {
        sp_sites <- NULL
    }

    message(paste("Fitting occupancy model to species", sp_code, "for year", year_sel, Sys.time()))

    # Fit model
    if(config$package == "spOccupancy"){
        fitSpOccu(site_data_year, visit_data_year, config, spatial, sp_sites, sp_code, year_sel)
    } else if(config$package == "occuR"){
        fitOccuR(site_data_year, visit_data_year, config, spatial, sp_sites, verbose = varargs$print_fitting)
    }


}
