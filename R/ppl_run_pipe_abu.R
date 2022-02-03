#' Run abundance indicators pipeline
#'
#' @param sp_code SAFRING code of the species to run the pipeline for
#' @param year Year to run to the pipeline for
#' @param site Code for site of interest.
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambJAGS}
#' @param steps A character vector containing the steps of the pipeline to run.
#' Can contain: "fit", "summ". Defaults to all of them.
#' @param ... Other arguments passed on to other functions
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_abu <- function(sp_code, year, site, config,
                             steps = c("fit", "summ"), ...){

    if("fit" %in% steps){
        ppl_fit_ssm_model(sp_code, site, year, config, ...)
    }

    if("summ" %in% steps){
        ppl_summarize_ssm(sp_code, site, year, config, ...)
    }

}
