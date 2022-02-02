#' Run distribution indicators pipeline
#'
#' @param sp_code SAFRING code of the species to run the pipeline for
#' @param sp_name Common name of the species to run the pipeline for. This is
#' necessary for plots and summaries.
#' @param year Year to run to the pipeline for
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambOccuR}
#' @param steps A character vector containing the steps of the pipeline to run.
#' Defaults to all of them.
#' @param ... Other arguments passed on to other functions
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_distr <- function(sp_code, sp_name, year, config,
                               steps = c("data", "fit", "summ", "indtr"), ...){

    if("data" %in% steps){
        ppl_create_site_visit(sp_code, year, config, ...)
    }

    if("fit" %in% steps){
        ppl_fit_occur_model(sp_code, year, config, ...)
    }

    if("summ" %in% steps){
        ppl_summarize_occur(sp_code, sp_name, year, config, ...)
    }

    if("indtr" %in% steps){
        # Estimate area of occupancy (AOO) for the year or also previous years
        # if year is earlier than 2010 (i.e. initial time window)
        for(t in seq_along(config$years)){

            if((t > config$dyear/2) | (config$year < (2009 + config$dyear/2))){
                ppl_estimate_aoo(sp_code, year = config$years[t], config, ...)
            }

        }

        # Estimate annual change in AOO and if year is greater than 2017,
        # also estimate short term change in AOO (i.e. in previous 10 years)
        if(year != 2017){
            for(t in seq_along(config$years[-1])){

                if(config$years[t+1] > 2008){
                    ppl_estimate_daoo(sp_code, year = config$years[t+1], config, term = "annual", ...)
                }

            }

            for(t in seq_along(config$years)){

                if(config$years[t] > 2017){
                    ppl_estimate_daoo(sp_code, year = config$years[t], config, term = "short", ...)
                }

            }
        }
    }

}
