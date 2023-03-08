#' Fit spOccupancy model
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
#' @param ... Other arguments that might be needed (e.g. for messages)
#'
#' @return Either a spOccupancy model fit or the integer 3, indicating that model fit
#' failed.
#' @export
#'
#' @examples
fitSpOccu <- function(site_data_year, visit_data_year, config, spatial = FALSE, sp_sites, ...){

    # Prepare data for spOccupancy
    occu_data <- prepSpOccuData_single(site_data_year, visit_data_year, config, spatial = spatial, sp_sites)


    # Define models -----------------------------------------------------------

    # Priors and initial values
    # Note: we could use posterior of previous years models to define priors

    # For single season models
    # Number of samples
    n_samples <- 2e4
    batch_length <- 25
    n_batch  <- n_samples/batch_length

    if(spatial){

        # Pair-wise distances between all sites
        dist_sites <- stats::dist(occu_data$coords)/1000

        # Specify list of inits
        inits <- list(alpha = 0,
                      beta = 0,
                      z = apply(occu_data$y, 1, max, na.rm = TRUE),
                      sigma.sq = 1,
                      phi = 3 / mean(dist_sites),
                      w = rep(0, nrow(occu_data$y)))

        # Priors
        priors <- list(alpha.normal = list(mean = 0, var = 2.72),
                       beta.normal = list(mean = 0, var = 2.72),
                       sigma.sq.ig = c(2, 1),
                       phi.unif = c(3/(1000*min(dist_sites)), 3/(min(dist_sites))))

        # Run model
        fit <- spOccupancy::spPGOcc(occ.formula = reformulate(config$site_mod),
                                    det.formula = reformulate(config$visit_mod),
                                    cov.model = "exponential", NNGP = TRUE, n.neighbors = 10,
                                    data = occu_data, inits = inits, priors = priors,
                                    batch.length = batch_length, n.batch = n_batch, n.burn = 2000,
                                    accept.rate = 0.43, tuning = list(phi = 4),
                                    n.omp.threads = 3, n.thin = 20, n.chains = 3,
                                    verbose = TRUE, n.report = 200)

    } else {

        filename <- paste0("occu_fit_", config$package, "_", year_sel-1, "_", sp_code, ".rds")

        if(file.exists(file.path(config$out_dir, sp_code, filename))){

            # Load previous fit
            prev_fit <- readRDS(file.path(config$out_dir, sp_code, filename))

            # Define priors
            priors <- defineSpOccupancyPriors(prev_fit)

            # Specify list of inits
            inits <- list(alpha = priors$alpha.normal$mean,
                          beta = priors$beta.normal$mean,
                          z = apply(occu_data$y, 1, max, na.rm = TRUE))

        } else {

            # Generic priors
            priors <- list(alpha.normal = list(mean = 0, var = 2.5),
                           beta.normal = list(mean = 0, var = 2.5))

            # Specify list of inits
            inits <- list(alpha = 0,
                          beta = 0,
                          z = apply(occu_data$y, 1, max, na.rm = TRUE))

        }


        # Run model

        fit <- tryCatch({
            out <- spOccupancy::PGOcc(occ.formula = reformulate(config$occ_mod),
                                      det.formula = reformulate(config$det_mod),
                                      data = occu_data, inits = inits, priors = priors,
                                      n.samples = n_samples, n.omp.threads = 1,
                                      n.thin = 20, n.chains = 3,
                                      verbose = TRUE, n.report = n_samples)

            out

        },
        error = function(cond) {
            sink(file.path(config$out_dir, paste0("reports/error_occu_fit_", year_sel, "_", sp_code, ".rds")))
            print(cond)
            sink()
            message(cond)
            return(NULL)
        },
        warning = function(cond) {
            sink(file.path(config$out_dir, paste0("reports/warning_occu_fit_", year_sel, "_", sp_code, ".rds")))
            print(cond)
            sink()
            message(cond)
            return(out)
        })

    }

    # Save fit and return 0 if success
    if(!is.null(fit)){

        # Save covariate scaling factors
        fit$det.scale <- list(scale = unlist(lapply(occu_data$det.covs, attr, "scaled:scale")),
                              center = unlist(lapply(occu_data$det.covs, attr, "scaled:center")))

        fit$occ.scale <- list(scale = attr(occu_data$occ.covs, "scaled:scale"),
                              center = attr(occu_data$occ.covs, "scaled:center"))


        return(fit)

    } else {

        return(3)

    }

}
