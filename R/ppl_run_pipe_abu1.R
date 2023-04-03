#' Run abundance pipeline ABU1
#'
#' @param sp_code SAFRING reference number of the species we want to analyze.
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambJAGS}
#' @param steps Pipeline steps to run. It can be one or more of: c("data", "fit", "diagnose", "summary").
#' @param prep_data_steps Data preparation steps to pass on to \link{ppl_create_data_ssm}
#' @param ... Other parameters to pass on to \link{prepGEESpCountData}
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_abu1 <- function(sp_code, config, steps = c("data", "fit", "diagnose", "summary"),
                              prep_data_steps, ...){

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
                 year = config$years_ch, data = NA, fit = NA, diagnose = NA, summary = NA,
                 package = config$package, notes = NA)


    # Prepare data ------------------------------------------------------------


    if("data" %in% steps){
        counts <- ppl_create_data_ssm(sp_code = sp_code, year = config$year,
                                      config = config,
                                      steps = prep_data_steps, ...)

        if(is.numeric(counts) && counts == 1){
            return(1)
        }

        # Save counts to disk
        outfile <- setSpOutFilePath("abu_model_data", config, sp_code, ".csv")
        utils::write.csv(counts, outfile, row.names = FALSE)
        message(paste("Final counts dataset saved at", outfile))

        # Log data status
        ppl_log["data"] <- 0
    }


    # Fit model ---------------------------------------------------------------

    if("fit" %in% steps){

        # Read data in if no counts are found in the environment
        countsfile <- setSpOutFilePath("abu_model_data", config, sp_code, ".csv")
        if(!exists("counts") && file.exists(countsfile)){
            counts <- utils::read.csv(countsfile)
        }

        if(exists("counts")){

            # Fit model
            fit <- ppl_fit_ssm_model(counts, sp_code, config)

            # Save
            saveRDS(fit, setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))

            # Log fit status
            ppl_log["fit"] <- 0

        } else {
            # Log fit status
            ppl_log["fit"] <- 1
        }
    }


    # Diagnose fit ------------------------------------------------------------

    if("diagnose" %in% steps){

        message(paste0("Diagnostics on species ", sp_code, " SSM for years ", paste(config$year_range, collapse = "-")))

        # Skip species if less than 5 suitable sites were detected during fitting model fitting
        error_file <- setSpOutFilePath("Less_5_sites", config, sp_code, ".txt")

        if(file.exists(error_file)){
            warning(paste0("Less than 5 sites for species ", sp_code, " ", paste(config$year_range, collapse = "-"), ". Cannot diagnose"))
        } else {

            # Proceed with diagnostics

            # Load model fit
            fit <- readRDS(setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))
            fit_stats <- BIRDIE:::processJAGSoutput(fit, DIC = FALSE, params.omit = NULL)
            rm(fit) # to free memory

            # Load count data
            counts <- read.csv(setSpOutFilePath("abu_model_data", config, sp_code, ".csv"))

            diags <- ppl_diagnose_ssm(fit_stats, counts, sp_code, config)

            diagfile <- setSpOutFilePath("abu_diagnostics", config, sp_code, ".csv")
            utils::write.csv(diags, diagfile, row.names = FALSE)

            # Log diagnosis status
            diag_log <- diags[,c("nc_pars", "Tmean", "Tsd", "Tdiff")]
            diag_log <- paste(diag_log, collapse = "-")
            ppl_log["diagnose"] <- diag_log

        }

    }

    # Summarise estimates -----------------------------------------------------

    if("summary" %in% steps){

        # Load counts and fit
        countsfile <- setSpOutFilePath("abu_model_data", config, sp_code, ".csv")
        if(!exists("counts") && file.exists(countsfile)){
            counts <- utils::read.csv(countsfile)
        }

        fitfile <- setSpOutFilePath("ssm_fit", config, sp_code, ".rds")
        if(!exists("fit") && file.exists(fitfile)){
            fit <- readRDS(fitfile)
        }

        # If files are available summarise fit
        if(exists("counts") && exists("fit")){

            summs <- ppl_summarise_ssm(fit, counts, sp_code, config)

            # Save counts to disk
            outfile <- setSpOutFilePath("ssm_pred", config, sp_code, "_all.csv")
            utils::write.csv(summs, outfile, row.names = FALSE)
            message(paste("Summaries saved at", outfile))

            # Log summary status
            ppl_log["summary"] <- 0

        } else {
            # Log summary status
            ppl_log["summary"] <- 1
        }

    }

    # More uninteresting activity log -----------------------------------------

    # Close log
    if(length(logfile) == 0) logfile <- NULL
    createLog(config, logfile, full_log = ppl_log)

    return(0)

}
