#' Run distribution indicators pipeline branch 1
#'
#' @description Branch 1 of the distribution indicators pipeline estimates
#' occupancy probabilities in South Africa for a selected species. These
#' occupancy probabilities form the basis for building more elaborated
#' indicators in other pipeline branches.
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

    # Find most recent log file
    logfile <- file.info(list.files(file.path(config$out_dir, "reports"), full.names = TRUE)) %>%
        dplyr::mutate(file = row.names(.)) %>%
        dplyr::filter(grepl("pipe_log", file)) %>%
        dplyr::arrange(desc(ctime)) %>%
        dplyr::slice(1) %>%
        dplyr::pull(file)

    # Open log with current time
    s <- Sys.time()
    attr(s,"tzone") <- "Africa/Johannesburg"

    ppl_log <- c(date_time = format(s), species = sp_code, model = "occ",
                 data = NA, fit = NA, diagnose = NA, summary = NA, notes = NA)

    if("data" %in% steps){
        ppl_create_site_visit(sp_code, config, ...)
        ppl_log["data"] <- 0
    }

    if("fit" %in% steps){

        fit_out <- ppl_fit_occu_model(sp_code, year_sel = year, config, ...)

        # Generate reports
        fit_status <- logFitStatus(fit_out, year, sp_code, config)

        # Create log
        ppl_log["fit"] <- fit_status
        createLog(config, logfile, full_log = ppl_log)

        # Save fit if the process was successful or return the status otherwise
        if(fit_status == 0){
            saveRDS(fit_out, file.path(config$out_dir, sp_code, paste0("occu_fit_", year, "_", sp_code, ".rds")))
        } else {
            return(fit_status)
        }

        # set pipeline status
        ppl_status <- fit_status

    }

    if("diagnose" %in% steps){

        fit <- readRDS(file.path(config$out_dir, sp_code, paste0("occu_fit_", year, "_", sp_code, ".rds")))

        diag_out <- diagnoseSpOccu(fit, sp_code, config, year)

        saveRDS(diag_out,
            file.path(config$out_dir, sp_code, paste0("occu_ppc_", year, "_", sp_code, ".rds")))

        # Create log
        ppl_log["diagnose"] <- mean(diag_out$fit.y.rep > diag_out$fit.y) # Bayes p

        rm(diag_out)

    }

    if("summary" %in% steps){

        if(!exists("fit")){
            fit <- readRDS(file.path(config$out_dir, sp_code, paste0("occu_fit_", year, "_", sp_code, ".rds")))
        }
        ppl_summarise_occu(fit, sp_code, sp_name, year, config)

        # Create log
        ppl_log["summary"] <- 0

    }

    # Close log
    createLog(config, logfile, full_log = ppl_log)

    if(exists("ppl_status") && ppl_status != 0){
        return(ppl_status)
    } else {
        return(0)
    }


}
