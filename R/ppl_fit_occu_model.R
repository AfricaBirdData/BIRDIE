#' Fit occupancy model
#'
#' @inheritParams ppl_run_pipe_dst1
#'
#' @return
#' @import occuR
#' @import mgcv
#' @export
#'
#' @examples
ppl_fit_occu_model <- function(sp_code, year, config, ...){

    varargs <- list(...)

    # Prepare data ------------------------------------------------------------

    # File names
    visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
    sitefile <- file.path(config$out_dir, paste0("occu_site_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

    # Read in site and visit data
    site_data <- utils::read.csv(sitefile)
    visit_data <- utils::read.csv(visitfile)
    det_data <- utils::read.csv(detfile)

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

    # Add detection info to visit data
    visit_data <- visit_data %>%
        dplyr::left_join(det_data,
                         by = c("CardNo", "StartDate", "year", "Pentad")) %>%
        dplyr::mutate(Spp = ifelse(obs == 0, "-", 1))

    for(t in seq_along(config$years)){

        year_sel <- config$years[t]

        # Subset data sets
        site_data_year <- site_data %>%
            dplyr::filter(year == year_sel)

        visit_data_year <- visit_data %>%
            dplyr::filter(year == year_sel)

        # Prepare for spOccupancy
        message(paste("Fitting occupancy model to species", sp_code, "for year", year_sel, Sys.time()))
        occu_data <- prepSpOccuData_single(site_data_year, visit_data_year, config)


        # Define models -----------------------------------------------------------

        spatial <- ifelse(t == 2, TRUE, FALSE)

        # Detection covariates
        visit_mod <- c("(1|obs_id)", "(1|site_id)", "log(hours+1)", "prcp", "tdiff", "cwac")

        # Occupancy covariates
        site_mod <- c("log_dist_coast", "watext", "log_watext", "watrec", "ndvi", "elev",
                      "prcp", "tdiff")

        # Priors and initial values
        # Note: we could use posterior of previous years models to define priors

        # For single season models
        # Number of samples
        n_samples <- 2e4
        batch_length <- 25
        n_batch  <- n_samples/batch_length

        if(spatial){

            # Add site coordinates to data
            pentads <- ABAP::getRegionPentads("country", "South Africa") %>%
                dplyr::filter(Name %in% unique(site_data_year$Name))

            sf::st_agr(pentads) = "constant"
            aeaproj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

            pentads <- sf::st_transform(pentads, aeaproj)

            occu_data$coords <- pentads %>%
                dplyr::arrange(pentad) %>%
                sf::st_centroid() %>%
                sf::st_coordinates()

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
            fit <- spPGOcc(occ.formula = reformulate(c(site_mod, "watrec*watext")),
                           det.formula = reformulate(visit_mod),
                           cov.model = "exponential", NNGP = TRUE, n.neighbors = 10,
                           data = occu_data, inits = inits, priors = priors,
                           batch.length = batch_length, n.batch = n_batch, n.burn = 2000,
                           accept.rate = 0.43, tuning = list(phi = 4),
                           n.omp.threads = 6, n.thin = 20, n.chains = 3,
                           verbose = TRUE, n.report = 200)

        } else {

            # Specify list of inits
            inits <- list(alpha = 0,
                          beta = 0,
                          z = apply(occu_data$y, 1, max, na.rm = TRUE))


            # Priors
            priors <- list(alpha.normal = list(mean = 0, var = 2.5),
                           beta.normal = list(mean = 0, var = 2.5))


            # Run model

            fit <- tryCatch({
                out <- spOccupancy::PGOcc(occ.formula = reformulate(c(site_mod, "watrec*watext")),
                                   det.formula = reformulate(visit_mod),
                                   data = occu_data, inits = inits, priors = priors,
                                   n.samples = n_samples, n.omp.threads = 6,
                                   n.thin = 20, n.chains = 3,
                                   verbose = TRUE, n.report = 5000)

                out

            },
            error = function(cond) {
                sink(file.path(config$out_dir, paste0("reports/error_occu_fit_", year_sel, "_", sp_code, ".rds")),
                     split = TRUE)
                print(cond)
                sink()
                return(NULL)
            },
            warning = function(cond) {
                sink(file.path(config$out_dir, paste0("reports/warning_occu_fit_", year_sel, "_", sp_code, ".rds")),
                     split = TRUE)
                print(cond)
                sink()
                return(out)
            })

        }


        # Save fit and return 0 if success
        if(!is.null(fit)){
            saveRDS(fit, file.path(config$out_dir, sp_code, paste0("occu_fit_", year_sel, "_", sp_code, ".rds")))
        }
    }

    return(0)

}
