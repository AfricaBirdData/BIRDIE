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
#' See \link{configPreambOccuR}
#' @param steps A character vector containing the steps of the pipeline to run.
#' Can contain: "data", "fit", "summ". Defaults to all of them.
#' @param ... Other arguments passed on to other functions
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_dst1 <- function(sp_code, sp_name, year, config,
                              steps = c("data", "fit", "summ"), ...){

    if("data" %in% steps){
        ppl_create_site_visit(sp_code, year, config, ...)
    }

    if("fit" %in% steps){

        fit_status <- ppl_fit_occur_model(sp_code, year, config, ...)

        # Stop if there are no detections
        if(fit_status == 1){
            warning(paste("There are no detections for species", sp_code))
            return(1)
        } else if(fit_status == 2){
            warning(paste("Species", sp_code, "detected in less than 5 pentads"))
            return(1)
        } else {
            # set pipeline status
            ppl_status <- fit_status
        }
    }

    if("summ" %in% steps){
        ppl_summarize_occur(sp_code, sp_name, year, config, ...)
    }

    return(ppl_status)

}
