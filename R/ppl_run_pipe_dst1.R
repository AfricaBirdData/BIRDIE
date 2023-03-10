#' Run distribution indicators pipeline module 1
#'
#' @description Module 1 of the distribution indicators pipeline estimates
#' occupancy probabilities in South Africa for a selected species. These
#' occupancy probabilities form the basis for building more elaborated
#' indicators in other pipeline modules.
#' @param sp_code SAFRING code of the species to run the pipeline for
#' @param sp_name Common name of the species to run the pipeline for. This is
#' necessary for plots and summaries.
#' @param year Year to run to the pipeline for
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambOccu}
#' @param steps A character vector containing the steps of the pipeline to run.
#' Can contain: "data", "fit", "diagnose", "summary". Defaults to all of them.
#' @param ... Other arguments passed on to other functions
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_dst1 <- function(sp_code, sp_name, year, config,
                              steps = c("data", "fit", "diagnose", "summary"), ...){

    varargs_dst1 <- (...)

    # Uninteresting activity logs ---------------------------------------------

    # Find most recent log file
    logfile <- file.info(list.files(file.path(config$out_dir, "reports"), full.names = TRUE)) %>%
        dplyr::mutate(file = row.names(.)) %>%
        dplyr::filter(grepl(paste0("pipe_log_", config$module), file)) %>%
        dplyr::arrange(desc(ctime)) %>%
        dplyr::slice(1) %>%
        dplyr::pull(file)

    # Open log with current time
    s <- Sys.time()
    attr(s,"tzone") <- "Africa/Johannesburg"

    ppl_log <- c(date_time = format(s), species = sp_code, model = "occ",
                 year = year, data = NA, fit = NA, diagnose = NA, summary = NA,
                 package = config$package, notes = NA)


    # Download and prepare data for model fitting -----------------------------

    if("data" %in% steps){
        ppl_create_site_visit(sp_code, config, ...)

        # Log data status
        ppl_log["data"] <- 0
    }


    # Fit occupancy model -----------------------------------------------------

    if("fit" %in% steps){

        fit_out <- ppl_fit_occu_model(sp_code, year_sel = year, config, varargs_dst1$spatial, ...)

        # Generate reports
        fit_status <- logFitStatus(fit_out, year, sp_code, config)

        # Log fit status
        ppl_log["fit"] <- fit_status

        # Save fit if the process was successful or return the status otherwise
        if(fit_status == 0){
            filename <- paste0("occu_fit_", config$package, "_", year, "_", sp_code, ".rds")
            saveRDS(fit_out, file.path(config$out_dir, sp_code, filename))
        } else {
            return(fit_status)
        }

        # set pipeline status
        ppl_status <- fit_status

    }


    # Diagnose model fit ------------------------------------------------------

    if("diagnose" %in% steps){

        filename <- paste0("occu_fit_", config$package, "_", year, "_", sp_code, ".rds")
        fit <- readRDS(file.path(config$out_dir, sp_code, filename))

        if(config$package == "spOccupancy"){
            diag_out <- diagnoseSpOccu(fit, sp_code, config, year)
        } else if(config$package == "occuR"){
            diag_out <- diagnoseOccuR(fit, sp_code, config, year)
        }

        saveRDS(diag_out,
                file.path(config$out_dir, sp_code, paste0("occu_ppc_", config$package, "_", year, "_", sp_code, ".rds")))

        # Create log
        ppl_log["diagnose"] <- mean(diag_out$fit.y.rep > diag_out$fit.y) # Bayes p

        rm(diag_out)

    }


    # Summarise predictions ---------------------------------------------------

    if("summary" %in% steps){

        if(!exists("fit")){
            filename <- paste0("occu_fit_", config$package, "_", year, "_", sp_code, ".rds")
            fit <- readRDS(file.path(config$out_dir, sp_code, filename))
        }

        ppl_summarise_occu(fit, sp_code, sp_name, year, config)

        # Create log
        ppl_log["summary"] <- 0

    }


    # More uninteresting activity log -----------------------------------------


    # Close log
    if(length(logfile) == 0) logfile <- NULL
    createLog(config, logfile, full_log = ppl_log)

    if(exists("ppl_status") && ppl_status != 0){
        return(ppl_status)
    } else {
        return(0)
    }


}
