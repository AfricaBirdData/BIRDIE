#' Run abundance pipeline ABU1
#'
#' @param sp_code SAFRING reference number of the species we want to analyze.
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambJAGS}
#' @param steps Pipeline steps to run. It can be one or more of: c("data", "fit", "summary").
#' @param prep_data_steps Data preparation steps to pass on to \link{ppl_create_data_ssm}
#' @param ... Other parameters to pass on to \link{prepGEESpCountData}
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_abu1 <- function(sp_code, config, steps = c("data", "fit", "summary"),
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
                 year = year, data = NA, fit = NA, diagnose = NA, summary = NA,
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
        utils::write.csv(counts_mod, outfile, row.names = FALSE)
        message(paste("Final counts dataset saved at", outfile))

        # Log data status
        ppl_log["data"] <- 0
    }


    # Fit model ---------------------------------------------------------------

    if("fit" %in% steps){

        # Read data in if no counts are found in the environment
        if(!exists(counts)){
            counts <- utils::read.csv(setSpOutFilePath("abu_model_data", config, sp_code, ".csv"))
        }

        # Fit model
        fit <- ppl_fit_ssm_model(counts, sp_code, config)

        # Save
        saveRDS(fit, setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))

        # Log fit status
        ppl_log["fit"] <- 0
    }


    if("summary" %in% steps){

        # Load counts and fit
        if(!exists(counts)){
            counts <- utils::read.csv(setSpOutFilePath("abu_model_data", config, sp_code, ".csv"))
        }
        if(!exists(fit)){
            fit <- readRDS(setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))
        }

        ppl_summarise_ssm(fit, counts, sp_code, config)

        # Log summary status
        ppl_log["summary"] <- 0

    }

    return(0)

}
