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

    occuRdata <- ppl_prep_occur_data(sp_code, year, config, ...)

    # Stop if there are no detections
    if(is.numeric(occuRdata) && occuRdata == 1){
        return(1)
    }

    # Model fitting -----------------------------------------------------------

    print(paste0("Fitting model at ", Sys.time(), ". This will take a while..."))

    visit_mod <- config$visit_mods

    # Determine whether the non-linear effect of month in p is necessary with the
    # simplest model for occupancy probabilities
    fit <- occuR::fit_occu(forms = list(reformulate(visit_mod, response = "p"),
                                        reformulate(config$site_mods[[3]], response = "psi")),
                           visit_data = occuRdata$visit,
                           site_data = occuRdata$site,
                           print = varargs$print_fitting)

    # Calculate degrees of freedom
    dof <- occuR::dof.occuR(fit, each = TRUE)

    # If non-linear effect of month on detection has very few dof,
    # fit linear effect to avoid singular covariance matrix
    print(paste(round(dof$p, 1), "DOF for the effect of month on p"))

    if(dof$p < 3){

        visit_mod <- c("1", "log(TotalHours+1)")

        # Create notification
        sink(file.path(config$fit_dir, sp_code, paste0("no_month_effect_", sp_code,".txt")))
        print("Model fitted without effect of month", split = TRUE)
        sink()
    }


    # Fit models sequentially if they don't work
    success <- FALSE
    m <- 0
    while(!success && m <= (length(config$site_mods) - 1)){

        m <- m + 1

        site_mod <- config$site_mods[[m]]
        print(paste("Trying model", m))

        tryCatch({

            fit <- occuR::fit_occu(forms = list(reformulate(visit_mod, response = "p"),
                                                reformulate(site_mod, response = "psi")),
                                   visit_data = occuRdata$visit,
                                   site_data = occuRdata$site,
                                   print = varargs$print_fitting)

            # Check if degrees of freedom can be calculated
            occuR::dof.occuR(fit, each = TRUE)

            success <- TRUE
            saveRDS(fit, file.path(config$fit_dir, sp_code, paste0("occur_fit_", config$years_ch, "_", sp_code, ".rds")))

        }, error = function(e){

            success <- FALSE

            sink(file.path(config$fit_dir, sp_code, paste0("failed_fit_", m, "_", sp_code,".txt")))
            print(e, split = TRUE)
            sink()

            }) # TryCatch fit
    }

    # Return 0 if success
    if(success){
        return(0)
    } else {
        return(2)
    }

}
