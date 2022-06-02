#' Fit occuR model
#'
#' @inheritParams ppl_run_pipe_distr
#'
#' @return
#' @import occuR
#' @import mgcv
#' @export
#'
#' @examples
ppl_fit_occur_model <- function(sp_code, year, config, ...){

    varargs <- list(...)

    # Prepare data ------------------------------------------------------------

    occuRdata <- prepOccuRData(sp_code, year, config, ...)

    # Stop if there are no detections
    if(is.numeric(occuRdata) && occuRdata == 1){
        return(1)
    }

    if(is.numeric(occuRdata) && occuRdata == 2){
        return(2)
    }

    # Define models -----------------------------------------------------------

    # Define spatio-temporal effect
    if(config$dur > 2){
        sptemp <- paste0("t2(lon, lat, occasion, k = c(", config$dim_grid, ", ", config$dur, "), bs = c('ts', 'cs'), d = c(2, 1))")
    } else {
        sptemp <- paste0("t2(lon, lat, k = ", config$dim_grid ,", bs = 'ts')")
    }

    # Detection
    visit_mod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

    # Occupancy
    site_mods <- list(mod1 = c("-1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec", sptemp),
                      mod2 = c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"),
                      mod3 = c("1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec"))


    # Model fitting -----------------------------------------------------------

    message(paste0("Fitting model at ", Sys.time(), ". This will take a while..."))

    # Determine whether the non-linear effect of month in p is necessary with the
    # simplest model for occupancy probabilities
    fit <- occuR::fit_occu(forms = list(reformulate(visit_mod, response = "p"),
                                        reformulate(site_mods[[3]], response = "psi")),
                           visit_data = occuRdata$visit,
                           site_data = occuRdata$site,
                           print = varargs$print_fitting)

    # Calculate degrees of freedom
    dof <- occuR::dof.occuR(fit, each = TRUE)

    # If non-linear effect of month on detection has very few dof,
    # fit linear effect to avoid singular covariance matrix
    message(paste(round(dof$p, 1), "DOF for the effect of month on p"))

    if(dof$p < 3){

        visit_mod <- c("1", "log(TotalHours+1)")

        # Create notification
        sink(file.path(config$out_dir, sp_code, paste0("no_month_effect_", sp_code,".txt")))
        message("Model fitted without effect of month")
        sink()
    }


    # Fit models sequentially if they don't work
    success <- FALSE
    m <- 0
    while(!success && m <= (length(site_mods) - 1)){

        m <- m + 1

        site_mod <- site_mods[[m]]
        message(paste("Trying model", m))

        tryCatch({

            fit <- occuR::fit_occu(forms = list(reformulate(visit_mod, response = "p"),
                                                reformulate(site_mod, response = "psi")),
                                   visit_data = occuRdata$visit,
                                   site_data = occuRdata$site,
                                   print = varargs$print_fitting)

            # Check if degrees of freedom can be calculated
            occuR::dof.occuR(fit, each = TRUE)

            success <- TRUE

        }, error = function(e){

            success <- FALSE

            sink(file.path(config$out_dir, sp_code, paste0("failed_fit_", m, "_", sp_code,".txt")))
            print(e, split = TRUE)
            sink()

            }) # TryCatch fit

        if(success && any(is.na(sqrt(diag(fit$res$cov.fixed))))){
            success <- FALSE
        }
    }

    # Save fit and return 0 if success
    if(success){
        saveRDS(fit, file.path(config$out_dir, sp_code, paste0("occur_fit_", config$years_ch, "_", sp_code, ".rds")))
        return(0)
    } else {
        return(2)
    }

}
