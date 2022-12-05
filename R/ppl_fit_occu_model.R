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

            # Pair-wise distances between all sites
            dist_sites <- dist(occu_data$coords)/1000

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
            out <- spPGOcc(occ.formula = reformulate(c(site_mod, "watrec*watext")),
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
            out <- spOccupancy::PGOcc(occ.formula = reformulate(c(site_mod, "watrec*watext")),
                                      det.formula = reformulate(visit_mod),
                                      data = occu_data, inits = inits, priors = priors,
                                      n.samples = n_samples, n.omp.threads = 10,
                                      n.thin = 20, n.chains = 3,
                                      verbose = TRUE, n.report = 5000)

        }


        # if(spatial){
        #
        #     # Pair-wise distances between all sites
        #     dist_sites <- dist(occu_data$coords)/1000
        #
        #     # Specify list of inits
        #     inits <- list(alpha = 0,
        #                   beta = 0,
        #                   z = apply(occu_data$y, 1, max, na.rm = TRUE),
        #                   sigma.sq = 1,
        #                   phi = 3 / mean(dist_sites),
        #                   w = rep(0, nrow(occu_data$y)))
        #
        #
        #     # Priors
        #     priors <- list(alpha.normal = list(mean = 0, var = 2.72),
        #                    beta.normal = list(mean = 0, var = 2.72),
        #                    sigma.sq.ig = c(2, 1),
        #                    phi.unif = c(3/(1000*min(dist_sites)), 3/(min(dist_sites))))
        #     # phi.unif = c(0.01, 1))
        #
        #     # source("functions/spPGOcc_tuned.R")
        #
        #     # Run model
        #     ptm <- proc.time()
        #     out <- spPGOcc(occ.formula = reformulate(c(site_mod, "watrec*watext")),
        #                    det.formula = reformulate(visit_mod),
        #                    cov.model = "exponential", NNGP = TRUE, n.neighbors = 10,
        #                    data = occu_data, inits = inits, priors = priors,
        #                    batch.length = 25, n.batch = 800, n.burn = 2000,
        #                    accept.rate = 0.43, tuning = list(phi = 4),
        #                    n.omp.threads = 6, n.thin = 20, n.chains = 3,
        #                    verbose = TRUE, n.report = 200)
        #     runtime <- proc.time() - ptm
        #
        # } else {
        #
        #     # Specify list of inits
        #     inits <- list(alpha = 0,
        #                   beta = 0,
        #                   sigma.sq.psi = 1, # occurrence random effect variances,
        #                   z = apply(occu_data$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0)))
        #
        #
        #     # Priors
        #     priors <- list(alpha.normal = list(mean = 0, var = 2.72),
        #                    beta.normal = list(mean = 0, var = 2.72),
        #                    sigma.sq.psi.ig = list(a = 0.1, b = 0.1))
        #
        #     # Number of samples
        #     n_samples <- 2e5
        #     batch_length <- 25
        #     n_batch  <- n_samples/batch_length
        #
        #     # Run model
        #     ptm <- proc.time()
        #     out <- tPGOcc(occ.formula = reformulate(c(site_mod, "watrec*watext")),
        #                   det.formula = reformulate(visit_mod),
        #                   ar1 = FALSE,
        #                   data = occu_data, inits = inits, priors = priors,
        #                   batch.length = batch_length, n.batch = n_batch, n.burn = 2000,
        #                   accept.rate = 0.43, tuning = list(phi = 4),
        #                   n.omp.threads = 6, n.thin = 20, n.chains = 3,
        #                   verbose = TRUE, n.report = 200)
        #     runtime <- proc.time() - ptm
        #
        # }

        # Diagnostics should go here. For now
        success <- TRUE

        # Save fit and return 0 if success
        if(success){
            saveRDS(fit, file.path(config$out_dir, sp_code, paste0("occu_fit_", year_sel, "_", sp_code, ".rds")))
            return(0)
        } else {
            return(3)
        }

    }
}






# OLD CODE ----------------------------------------------------------------







#     # Define models -----------------------------------------------------------
#
#     # Define spatio-temporal effect
#     if(config$dur > 2){
#         sptemp <- paste0("t2(lon, lat, occasion, k = c(", config$dim_grid, ", ", config$dur, "), bs = c('ts', 'cs'), d = c(2, 1))")
#     } else {
#         sptemp <- paste0("t2(lon, lat, k = ", config$dim_grid ,", bs = 'ts')")
#     }
#
#     # Detection
#     visit_mod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")
#
#     # Occupancy
#     site_mods <- list(mod1 = c("-1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec", sptemp),
#                       mod2 = c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"),
#                       mod3 = c("1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec"))
#
#
#     # Model fitting -----------------------------------------------------------
#
#     message(paste0("Fitting model at ", Sys.time(), ". This will take a while..."))
#
#     # Determine whether the non-linear effect of month in p is necessary with the
#     # simplest model for occupancy probabilities
#     fit <- occuR::fit_occu(forms = list(reformulate(visit_mod, response = "p"),
#                                         reformulate(site_mods[[3]], response = "psi")),
#                            visit_data = occu_data$visit,
#                            site_data = occu_data$site,
#                            print = varargs$print_fitting)
#
#     # Calculate degrees of freedom
#     dof <- occuR::dof.occuR(fit, each = TRUE)
#
#     # If non-linear effect of month on detection has very few dof,
#     # fit linear effect to avoid singular covariance matrix
#     message(paste(round(dof$p, 1), "DOF for the effect of month on p"))
#
#     if(dof$p < 3){
#
#         visit_mod <- c("1", "log(TotalHours+1)")
#
#         # Create notification
#         sink(file.path(config$out_dir, sp_code, paste0("no_month_effect_", sp_code, "_", year, ".txt")))
#         message("Model fitted without effect of month")
#         sink()
#     }
#
#
#     # Fit models sequentially if they don't work
#     success <- FALSE
#     m <- 0
#     while(!success && m <= (length(site_mods) - 1)){
#
#         m <- m + 1
#
#         site_mod <- site_mods[[m]]
#         message(paste("Trying model", m))
#
#         tryCatch({
#
#             fit <- occuR::fit_occu(forms = list(reformulate(visit_mod, response = "p"),
#                                                 reformulate(site_mod, response = "psi")),
#                                    visit_data = occu_data$visit,
#                                    site_data = occu_data$site,
#                                    print = varargs$print_fitting)
#
#             # Check if degrees of freedom can be calculated
#             occuR::dof.occuR(fit, each = TRUE)
#
#             success <- TRUE
#
#         }, error = function(e){
#
#             success <- FALSE
#
#             sink(file.path(config$out_dir, sp_code, paste0("failed_fit_", m, "_", sp_code,"_", year, ".txt")))
#             print(e, split = TRUE)
#             sink()
#
#             }) # TryCatch fit
#
#         if(success && any(is.na(sqrt(diag(fit$res$cov.fixed))))){
#             success <- FALSE
#         }
#     }
#
#     # Save fit and return 0 if success
#     if(success){
#         saveRDS(fit, file.path(config$out_dir, sp_code, paste0("occur_fit_", config$years_ch, "_", sp_code, ".rds")))
#         return(0)
#     } else {
#         return(2)
#     }
#
# }
