#' Run distribution indicators pipeline module 1
#'
#' @description Module 1 of the distribution indicators pipeline estimates
#' occupancy probabilities in South Africa for a selected species. These
#' occupancy probabilities form the basis for building more elaborated
#' indicators in other pipeline modules.
#' @param sp_code SAFRING code of the species to run the pipeline for
#' @param year Year to run to the pipeline for
#' @param config A list with pipeline configuration parameters.
#' See \code{\link{configPipeline}}
#' @param steps A character vector containing the steps of the pipeline to run.
#' Can contain: "data", "fit", "diagnose", "summary". Defaults to all of them.
#' @param ... Other arguments passed on to other functions
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_dst1 <- function(sp_code, year, config,
                              steps = c("data", "fit", "diagnose", "summary"), ...){

    varargs_dst1 <- list(...)

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

        occu_data <- ppl_create_site_visit(config, sp_code,
                                           force_gee_dwld = varargs_dst1$force_gee_dwld ,
                                           force_site_visit = varargs_dst1$force_site_visit,
                                           force_abap_dwld = varargs_dst1$force_abap_dwld,
                                           monitor_gee = varargs_dst1$monitor_gee)

        # Save data files
        visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
        sitefile <- file.path(config$out_dir, paste0("occu_site_dat_sa_", config$years_ch, ".csv"))
        detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

        # Save data
        occu_data$visit %>%
            dplyr::select(-obs) %>%
            utils::write.csv(visitfile, row.names = FALSE)

        occu_data$site %>%
            utils::write.csv(sitefile, row.names = FALSE)

        occu_data$visit %>%
            dplyr::select(CardNo, StartDate, Pentad, year, obs) %>%
            utils::write.csv(detfile, row.names = FALSE)

        message(paste("Site-visit data saved at", config$out_dir))

        # Preliminary data status
        data_status <- 0

        # Check if there are problems with the data and create logs
        if(!1 %in% unique(occu_data$visit$obs)){
            warning(paste("No detection of species", sp_code))
            # Log data status
            data_status <- 1
        }

        # Or species detected in too few Pentads
        min_pentads <- 2*length(config$occ_mod)

        year_sel <- year
        n_pentads <- occu_data$visit %>%
            dplyr::count(year, Pentad, obs) %>%
            dplyr::filter(year == year_sel, obs == 1) %>%
            nrow()
        rm(year_sel)

        if(n_pentads < min_pentads){
            warning(paste("Species", sp_code, "detected in less than", min_pentads, "pentads"))
            # Log data status
            data_status <- 2
        }

        # Log data status
        ppl_log["data"] <- data_status

        # Display and save logs
        if(data_status != 0){
            txt <- dplyr::case_when(data_status == 1 ~ paste0("no_detections_", year,"_", sp_code, ".txt"),
                                    data_status == 2 ~ paste0("less_than_", min_pentads, "_pentads_", year,"_", sp_code, ".txt"))

            conv_file <- file.path(config$out_dir, sp_code, txt)
            sink(conv_file)
            print(gsub(".txt", "", txt))
            sink()
            message(gsub(".txt", "", txt)) # to console

            createLog(config, logfile, full_log = ppl_log)

            return(data_status)
        }

    }


    # Fit occupancy model -----------------------------------------------------

    if("fit" %in% steps){

        # Dump a message for debugging
        sink(file.path(config$out_dir, "reports", paste0(Sys.time(), "_Fitting_species_", sp_code, "_year_", year, ".txt")))
        print(paste(Sys.time(), "Fitting occupancy model to species", sp_code, "for year", year_sel))
        sink()

        # Fit model
        fit_out <- ppl_fit_occu_model(sp_code, year_sel = year, config, varargs_dst1$spatial, ...)

        # Generate reports
        # fit_status <- logFitStatus(fit_out, year, sp_code, config)

        if(is.numeric(fit_out)){ # this immediately flags a problem
            fit_status <- fit_out
        } else {
            fit_status <- 0
        }

        # Log fit status
        ppl_log["fit"] <- fit_status

        # Save fit if the process was successful or return the status otherwise
        if(fit_status == 0){
            filename <- setSpOutFilePath("occu_fit", config, year, sp_code, ".rds")
            saveRDS(fit_out, filename)
        } else {
            createLog(config, logfile, full_log = ppl_log)
            return(fit_status)
        }

        # set pipeline status
        ppl_status <- fit_status

    }


    # Diagnose model fit ------------------------------------------------------

    if("diagnose" %in% steps){

        message(paste0("Diagnostics on species ", sp_code, " occupancy model for years ", year))

        # Skip if there were errors during fitting model fitting
        min_pentads <- 2*length(config$occ_mod)

        txt <- c(paste0("no_detections_", year,"_", sp_code, ".txt"),
                 paste0("less_than_", min_pentads, "_pentads_", year,"_", sp_code, ".txt"),
                 paste0("model_fit_failed_", year,"_", sp_code, ".txt"))

        error_files <- file.path(config$out_dir, sp_code, txt)

        if(any(file.exists(error_files))){
            warning(paste0("There were errors during model fitting for species ", sp_code, " ", year, ". Cannot diagnose"))
        } else {

            # Proceed with diagnostics

            # Load model fit
            filename <- setSpOutFilePath("occu_fit", config, year, sp_code, ".rds")
            fit <- readRDS(filename)

            diags <- ppl_diagnose_occu(fit, data=NULL, sp_code, year)

            # Save basic diagnostics (convergence and Bayesian p-value)
            diags$rhat$bayes_p <- diags$gof$bayes_p

            diagfile <- setSpOutFilePath("occu_diagnostics", config, year, sp_code, ".csv")
            utils::write.csv(diags$rhat, diagfile, row.names = FALSE)

            # Save posterior predictive checks
            ppcfile <- setSpOutFilePath("occu_ppc", config, year, sp_code, ".rds")
            saveRDS(diags$gof, ppcfile)

            # Log diagnosis status
            diag_log <- diags$rhat[,c("nc_pars", "bayes_p")]
            diag_log <- paste(diag_log, collapse = "-")
            ppl_log["diagnose"] <- diag_log

        }

    }


    # Summarise predictions ---------------------------------------------------

    if("summary" %in% steps){

        filename <- setSpOutFilePath("occu_fit", config, year, sp_code, ".rds")

        if(exists("fit")){
            fit <- fit
        } else if(file.exists(filename)){
            fit <- readRDS(filename)
        }

        if(exists("fit")){

            ppl_summarise_occu(fit, sp_code, year, config)
            # Create log
            ppl_log["summary"] <- 0

        } else {

            # If no model fit is found, then we use raw data to create the
            # predictions data frame.
            pred_from_data <- createPredFromAbap(sp_code, year, config)

            # Save predictions data frame
            summfile <- setSpOutFilePath("occu_pred", config, year, sp_code, ".csv")
            pred_from_data[["df"]] %>%
                write.csv(file.path(config$out_dir, sp_code, summfile),
                          row.names = FALSE)

            # Save plot
            summplot <- setSpOutFilePath("occu_occu", config, year, sp_code, ".png")
            ggplot2::ggsave(summplot, pred_from_data[["occu_plot"]])

            ppl_log["summary"] <- 1
        }

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
