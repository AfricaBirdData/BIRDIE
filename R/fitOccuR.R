#' Fit occuR occupancy model
#'
#' @param site_data_year Occupancy site data for a single year and species
#' (see \code{\link{ppl_create_site_visit}})
#' @param visit_data_year Occupancy visit data for a single year and species
#' (see \code{\link{ppl_create_site_visit}})
#' @param config A list with pipeline configuration parameters
#' (see \code{\link{configPreambOccu}}).
#' @param spatial Logical, indicating whether spatial random effects should be
#' included in the model (TRUE) or not (FALSE, default)
#' @param sp_sites Spatial object containing the pentads in `site_data_year`.
#' @param verbose Logical indicating whether information should be printed
#' during model fitting.
#'
#' @return Either an occuR model fit or the integer 3, indicating that model fit
#' failed.
#' @export
#'
#' @examples
fitOccuR <- function(site_data_year, visit_data_year, config, spatial = FALSE, sp_sites, verbose){

    # Prepare data for occuR
    occu_data <- prepOccuRData(site_data_year, visit_data_year, config, spatial = spatial, sp_sites)


    # Define models -----------------------------------------------------------

    # forms, visit_data, site_data, start = NULL, print = TRUE
    occu_data$visit <- occu_data$visit %>%
        dplyr::mutate(obs_id = factor(obs_id),
                      site_id = factor(site_id))

    # Run model
    fit <- tryCatch({
        occuR::fit_occu(forms = list(stats::reformulate(config$occ_mod, response = "psi"),
                                     stats::reformulate(config$det_mod, response = "p")),
                        visit_data = occu_data$visit,
                        site_data = occu_data$site,
                        print = verbose)

    },
    error = function(cond) {
        filename <- paste0("reports/error_occu_fit_", config$package, "_", year_sel, "_", sp_code, ".txt")
        sink(file.path(config$out_dir, filename))
        print(cond)
        sink()
        message(cond)
    })

    # Save fit and return 0 if success
    if(!is.null(fit)){

        # Save covariate scaling factors
        fit$det.scale <- list(scale = unlist(lapply(occu_data$visit, attr, "scaled:scale")),
                              center = unlist(lapply(occu_data$visit, attr, "scaled:center")))

        fit$occ.scale <- list(scale = unlist(lapply(occu_data$site, attr, "scaled:scale")),
                              center = unlist(lapply(occu_data$site, attr, "scaled:center")))


        return(fit)

    } else {

        return(3)

    }

}
