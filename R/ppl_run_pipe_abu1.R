#' Run abundance pipeline ABU1
#'
#' @param sp_code SAFRING reference number of the species we want to analyse.
#' @param config A list with pipeline configuration parameters.
#' See \code{\link{configPipeline}}
#' @param steps Pipeline steps to run. It can be one or more of: c("data", "fit", "diagnose", "summary").
#' @param prep_data_steps Data preparation steps to pass on to \code{\link{ppl_create_data_ssm}}
#' @param summary_scale Either "linear", summaries are given the linear scale or "model",
#' summaries are given in the modelling scale (typically log scale)
#' @param ... Other parameters to pass on to \code{\link{prepGEECatchmData}}
#'
#' @return This function will run the whole abundance module of the BIRDIE pipeline
#' @export
#'
#' @examples
ppl_run_pipe_abu1 <- function(sp_code, config, steps = c("data", "fit", "diagnose", "summary"),
                              prep_data_steps, summary_scale = c("linear", "model"), ...){

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

    ppl_log <- c(date_time = format(s), species = sp_code, model = "abu",
                 year = config$years_ch, data = NA, fit = NA, diagnose = NA, summary = NA,
                 package = config$package, notes = NA)


    # Prepare data ------------------------------------------------------------


    if("data" %in% steps){

        message(paste0(Sys.time()," Preparing SSM data for species ", sp_code, " for years ", paste(config$year_range, collapse = "-")))

        counts <- ppl_create_data_ssm(sp_code = sp_code, year = config$year,
                                      config = config,
                                      steps = prep_data_steps, ...)

        # if(is.numeric(counts) && counts == 1){
        #     return(1)
        # }

        # Save counts to disk
        if(!is.numeric(counts)){

            outfile <- setSpOutFilePath("abu_model_data", config, config$years_ch, sp_code, ".csv")
            utils::write.csv(counts, outfile, row.names = FALSE)
            message(paste("Final counts dataset saved at", outfile))

            # Log data status
            ppl_log["data"] <- 0

        } else {
            ppl_log["data"] <- counts
        }

    }


    # Fit model ---------------------------------------------------------------

    if("fit" %in% steps){

        message(paste0(Sys.time()," Fitting SSM for species ", sp_code, " for years ", paste(config$year_range, collapse = "-")))

        # Read data in if no counts are found in the environment
        countsfile <- setSpOutFilePath("abu_model_data", config, config$years_ch, sp_code, ".csv")

        # Remove any counts object to make sure we use the one on disk
        counts <- NULL

        if(file.exists(countsfile)){
            counts <- utils::read.csv(countsfile)
        }

        if(!is.null(counts)){

            # Fit model
            fit <- ppl_fit_ssm_model(counts, sp_code, config)

            # Save
            saveRDS(fit, setSpOutFilePath("ssm_fit", config, config$years_ch, sp_code, ".rds"))

            # Log fit status
            ppl_log["fit"] <- 0

        } else {
            warning("No data suitable for modelling found")
            # Log fit status
            ppl_log["fit"] <- 1
        }

    }


    # Diagnose fit ------------------------------------------------------------

    if("diagnose" %in% steps){

        message(paste0(Sys.time()," Diagnostics on species ", sp_code, " SSM for years ", paste(config$year_range, collapse = "-")))

        # Skip species if no model fit is found
        fit_file <- setSpOutFilePath("ssm_fit", config, config$years_ch, sp_code, ".rds")

        if(!file.exists(fit_file)){
            warning(paste0("No model fit for species ", sp_code, " ", paste(config$year_range, collapse = "-"), ". Cannot diagnose"))
        } else {

            # Proceed with diagnostics

            # Load model fit
            fit <- readRDS(fit_file)
            fit_stats <- BIRDIE::processJAGSoutput(fit, DIC = FALSE, params.omit = NULL)
            rm(fit) # to free memory

            # Load count data
            counts <- read.csv(setSpOutFilePath("abu_model_data", config, config$years_ch, sp_code, ".csv"))

            diags <- ppl_diagnose_ssm(fit_stats, counts, sp_code, config)

            diagfile <- setSpOutFilePath("abu_diagnostics", config, config$years_ch, sp_code, ".csv")
            utils::write.csv(diags, diagfile, row.names = FALSE)

            # Log diagnosis status
            diag_log <- diags[,c("nc_pars", "Tmean", "Tsd", "Tdiff")]
            diag_log <- paste(diag_log, collapse = "-")
            ppl_log["diagnose"] <- diag_log

        }

    }

    # Summarise estimates -----------------------------------------------------

    if("summary" %in% steps){

        message(paste0(Sys.time()," Summarising SSM data for species ", sp_code, " for years ", paste(config$year_range, collapse = "-")))

        # Load counts and fit
        countsfile <- setSpOutFilePath("abu_model_data", config, config$years_ch, sp_code, ".csv")

        if(file.exists(countsfile)){
            counts <- utils::read.csv(countsfile, colClasses = c(LocationCode = "character"))
        }

        fitfile <- setSpOutFilePath("ssm_fit", config, config$years_ch, sp_code, ".rds")
        if(file.exists(fitfile)){
            fit <- readRDS(fitfile)
        }

        # Extra counts at sites with insufficient data
        extracountsfile <- setSpOutFilePath("cwac_data_w_miss", config, config$years_ch, sp_code, ".csv")
        if(file.exists(extracountsfile)){
            extra_counts <- utils::read.csv(extracountsfile, colClasses = c(LocationCode = "character"))
        }

        # If files are available summarise fit
        if(exists("counts") && exists("fit")){

            summary_scale <- match.arg(summary_scale)
            linear <- TRUE
            if(summary_scale != "linear"){
                linear = FALSE
            }

            summs <- ppl_summarise_ssm(fit, counts, sp_code, linear=linear, config)

            # Add data from sites with insufficient data
            extra_counts <- read.csv(setSpOutFilePath("cwac_data_w_miss", config, config$years_ch, sp_code, ".csv"),
                                     colClasses = c(LocationCode = "character"))
            summs <- addExtraSitesToSummary(extra_counts, summs)

            # Save predictions to disk
            outfile <- setSpOutFilePath("ssm_pred", config, config$years_ch, sp_code, "_all.csv")
            utils::write.csv(summs, outfile, row.names = FALSE)
            message(paste("Summaries saved at", outfile))

            # Log summary status
            ppl_log["summary"] <- 0

            # If count and fit are not available return only counts

        } else if(exists("extra_counts")){

            summs <- extra_counts %>%
                dplyr::filter(Season %in% c("S", "W"))

            # But if there are not counts in summer or winter, then return empty
            if(nrow(summs) > 0){
                summs <- summs %>%
                    tidyr::pivot_wider(names_from = "Season", values_from = "count") %>%
                    dplyr::mutate(species = sp_code) %>%
                    dplyr::rename(site = "LocationCode",
                                  summer.count = "S",
                                  winter.count = "W") %>%
                    dplyr::select("species", "site", "year", "summer.count", "winter.count")
            } else {
                summs <- data.frame(species = sp_code,
                                    site = NA,
                                    year = config$years,
                                    summer.count = NA,
                                    winter.count = NA)
            }

            # Save predictions to disk
            outfile <- setSpOutFilePath("ssm_pred", config, config$years_ch, sp_code, "_all.csv")
            utils::write.csv(summs, outfile, row.names = FALSE)
            message(paste("Summaries saved at", outfile))

            # Log summary status
            ppl_log["summary"] <- 1
        } else {
            ppl_log["summary"] <- 1
        }

    }

    # More uninteresting activity log -----------------------------------------

    # Close log
    if(length(logfile) == 0) logfile <- NULL
    createLog(config, logfile, full_log = ppl_log)

    return(0)

}
