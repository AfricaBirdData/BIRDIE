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

    if("data" %in% steps){
        counts <- ppl_create_data_ssm(sp_code, config$year, config,
                                      steps = prep_data_steps, ...)

        if(is.numeric(counts) & counts == 1){
            return(1)
        }
    }


    # Fit model ---------------------------------------------------------------

    if("fit" %in% steps){
        ppl_fit_ssm_model(sp_code, config)
    }


    if("summary" %in% steps){
        ppl_summarize_ssm(sp_code, config)
    }

    return(0)

}
