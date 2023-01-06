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

    if("data" %in% steps){
        ppl_create_site_visit(sp_code, config, ...)
    }

    if("fit" %in% steps){

        fit_status <- ppl_fit_occu_model(sp_code, year_sel = year, config, ...)

        # Generate reports
        if(fit_status == 1){
            conv_file <- file.path(config$out_dir, reports, paste0("no_detections_", year,"_", sp_code, ".txt"))
            sink(conv_file)
            print(paste("no detections", year, sp_code))
            sink()
            message(paste("no detections", year, sp_code)) # to console
            return(fit_status)
        } else if(fit_status == 2){
            conv_file <- file.path(config$out_dir, reports, paste0("less_than_5_pentads_", year,"_", sp_code, ".txt"))
            sink(conv_file)
            print(paste("less than 5 pentads", year, sp_code))
            sink()
            message(paste("less than 5 pentads", year, sp_code)) # to console
            return(fit_status)
        } else if(fit_status == 3){
            conv_file <- file.path(config$out_dir, reports, paste0("model_fit_failed_", year,"_", sp_code, ".txt"))
            sink(conv_file)
            print(paste("model fit failed", year, sp_code))
            sink()
            message(paste("model fit failed", year, sp_code)) # to console
            return(fit_status)
        } else {
            saveRDS(fit_status, file.path(config$out_dir, sp_code, paste0("occu_fit_", year, "_", sp_code, ".rds")))
            fit_status <- 0
        }

        # set pipeline status
        ppl_status <- fit_status

    }

    if("diagnose" %in% steps){

        fit <- readRDS(file.path(config$out_dir, sp_code, paste0("occu_fit_", year, "_", sp_code, ".rds")))

        saveRDS(
            diagnoseSpOccu(fit, sp_code, config, year),
            file.path(config$out_dir, sp_code, paste0("occu_ppc_", year, "_", sp_code, ".rds"))
        )
    }

    if("summary" %in% steps){

        if(!exists("fit")){
            fit <- readRDS(file.path(config$out_dir, sp_code, paste0("occu_fit_", year, "_", sp_code, ".rds")))
        }
        ppl_summarise_occu(fit, sp_code, sp_name, year, config)

    }

    if(exists("ppl_status") && ppl_status != 0){
        return(ppl_status)
    } else {
        return(0)
    }


}
