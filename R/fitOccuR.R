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
    occu_data <- prepOccuRData(site_data_year, visit_data_year, config, spatial = spatial,
                               sp_sites, scale = TRUE, keep_sites = FALSE)


    # Define models -----------------------------------------------------------

    # forms, visit_data, site_data, start = NULL, print = TRUE
    occu_data$site <- occu_data$site %>%
        dplyr::mutate(site_id = factor(site_id))
    occu_data$visit <- occu_data$visit %>%
        dplyr::mutate(site_id = factor(site_id),
                      obs_id = factor(obs_id))

    # Run model
    tryCatch({
        fit <- occuR::fit_occu(forms = list(stats::reformulate(config$occ_mod, response = "psi"),
                                     stats::reformulate(config$det_mod, response = "p")),
                        visit_data = occu_data$visit,
                        site_data = occu_data$site,
                        print = verbose)

        if(fit$fit$convergence != 0){
            warning(paste(config$package, "model failed once to converge for species", sp_code, "and year", year_sel))
            fp <- fit$res$par.fixed
            rp <- fit$res$par.random
            start_vals <- list(beta_psi = fp[names(fp) == "beta_psi"],
                               beta_p = fp[names(fp) == "beta_p"],
                               z_psi = rp[names(rp) == "z_psi"],
                               z_p = rp[names(rp) == "z_p"],
                               log_lambda_psi = fp[names(fp) == "log_lambda_psi"],
                               log_lambda_p = fp[names(fp) == "log_lambda_p"],
                               gamma_psi = rp[names(rp) == "gamma_psi"],
                               gamma_p = rp[names(rp) == "gamma_p"],
                               lsig_gamma_psi = fp[names(fp) == "lsig_gamma_psi"],
                               lsig_gamma_p = fp[names(fp) == "lsig_gamma_p"])

            start_vals[sapply(start_vals, function(x) length(x) == 0)] <- 0

            fit <- occuR::fit_occu(forms = list(stats::reformulate(config$occ_mod, response = "psi"),
                                                stats::reformulate(config$det_mod, response = "p")),
                                   visit_data = occu_data$visit,
                                   site_data = occu_data$site,
                                   start = start_vals,
                                   print = verbose)

        }


        # if(fit$fit$convergence != 0){
        #     warning(paste(config$package, "model failed twice to converge for species", sp_code, "and year", year_sel))
        #     fit <- occuR::fit_occu(forms = list(stats::reformulate(config$occ_mod, response = "psi"),
        #                                         stats::reformulate(config$det_mod[-2], response = "p")),
        #                            visit_data = occu_data$visit,
        #                            site_data = occu_data$site,
        #                            print = verbose)
        #
        # }

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
