#' Fit occuR model
#'
#' @inheritParams ppl_run_pipe_distr
#'
#' @return
#' @import occuR
#' @export
#'
#' @examples
ppl_fit_occur_model <- function(sp_code, year, config, ...){

    varargs <- list(...)

    # Prepare data ------------------------------------------------------------

    occuRdata <- ppl_prep_occur_data(sp_code, year, config, ...)


    # Define models -----------------------------------------------------------

    # Detection
    visit_mod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

    # Occupancy
    site_mods <- list(mod1 = c("-1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec", config$sptemp),
                      mod2 = c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"),
                      mod3 = c("1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec"))


    # Model fitting -----------------------------------------------------------

    print(paste0("Fitting model at ", Sys.time(), ". This will take a while..."))

    # Fit models sequentially if they don't work
    success <- FALSE
    m <- 0
    while(!success && m <= length(site_mods)){
        m <- m + 1
        site_mod <- site_mods[[m]]
        print(paste("Trying model", m))
        tryCatch({
            fit <- occuR::fit_occu(forms = list(reformulate(visit_mod, response = "p"),
                                                reformulate(site_mod, response = "psi")),
                                   visit_data = occuRdata$visit,
                                   site_data = occuRdata$site,
                                   print = varargs$print_fitting)

            # Check if degrees of freedom can be calculated
            occuR::dof.occuR(fit)

            success <- TRUE
            saveRDS(fit, file.path(config$fit_dir, sp_code, paste0("occur_fit_", config$years_ch, "_", sp_code, ".rds")))

        }, error = function(e){
            success <- FALSE
            sink(file.path(config$fit_dir, sp_code, paste0("failed_fit_", m, "_", sp_code,".txt")))
            print(e)
            sink()
            }) # TryCatch fit
    }


}
