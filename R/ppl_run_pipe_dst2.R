#' Run distribution indicators pipeline branch 2
#'
#' @param sp_code SAFRING code of the species to run the pipeline for
#' @param indtr A character vector containing the indicators of the pipeline
#' to compute. Can contain: "aoo", "daoo". Defaults to all of them.
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambOccuR}
#' @param overwrite_indtr Logical. If TRUE, existing files in directories
#' corresponding to the species in config$species will be overwritten.
#' @param verbose Logical. If TRUE, the new line added to the indicator table
#' is displayed in the console.
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_dst2 <- function(sp_code, indtr = c("aoo", "daoo"), config,
                              overwrite_indtr, verbose, ...){

    # Create indicator file if it doesn't exist
    ppl_create_indtr_file(sp_code, config$year,
                          overwrite_indtr = overwrite_indtr)

    if("aoo" %in% indtr){
        # Estimate area of occupancy (AOO)
        aoo_status <- ppl_estimate_aoo(sp_code, config, verbose, ...)

        if(aoo_status == 1){
            print(paste("Too few detections to estimate AOO for species", sp_code))
            return(1)
        }

    }

    if("daoo" %in% indtr){
        # Estimate change in AOO
        ppl_estimate_daoo(sp_code, config, term = "annual", verbose,
                          overwrite = overwrite_indtr, ...)
    }

    return(0)

}
