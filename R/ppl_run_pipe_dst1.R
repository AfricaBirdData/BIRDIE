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

        # Stop if there are no detections
        if(fit_status == 1){
            sink(file.path(config$out_dir, sp_code,paste0("no_detections_", config$years_ch,"_", sp_code, ".txt")))
            sink()
            return(1)
        } else if(fit_status == 2){
            sink(file.path(config$out_dir, sp_code,paste0("less_than_5_pentads_", config$years_ch,"_", sp_code, ".txt")))
            sink()
            return(1)
        } else if(fit_status == 3){
            sink(file.path(config$out_dir, sp_code,paste0("model_fit_failed_", config$years_ch,"_", sp_code, ".txt")))
            sink()
            return(1)
        } else {
            # set pipeline status
            ppl_status <- fit_status
        }
    }

    if("diagnose" %in% steps){
        diagnoseSpOccu()
    }

    if("summary" %in% steps){
        ppl_summarize_occur(sp_code, sp_name, year, config, ...)
    }

    if(exists("ppl_status") && ppl_status != 0){
        return(ppl_status)
    } else {
        return(0)
    }


}
