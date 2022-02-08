#' Estimate distribution indicators for species type
#'
#' @param sp_type A character string with the species type for which
#' @param year Year for which indicators are to be estimated
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambOccuR}
#' @param force_predict If TRUE new prediction will be run from models and they
#' will be stored in a temporary directory. If FALSE cached predictions will be
#' used if they exist.
#' @param verbose Logical. If TRUE, the new line added to the indicator table
#' is displayed in the console.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ppl_estimate_distr_sp_type <- function(sp_type, year, config, force_predict = FALSE, verbose, ...){

    print("Estimating occupancy for group")

    estimateOccuSpType(sp_type, year, config, ...)

    print("Estimating Area Of Occupancy for group")

    for(t in seq_along(config$years)){
        estimateAOOSpType(sp_type, config$years[t], config, force_predict, verbose, ...)
    }

    print("Estimating change in Area Of Occupancy for group")

    # Estimate annual change in AOO and if year is greater than 2017,
    # also estimate short term change in AOO (i.e. in previous 10 years)
    for(t in seq_along(config$years[-1])){

        if(config$years[t+1] > 2008){
            ppl_estimate_daoo(sp_type, year = config$years[t+1], config, term = "annual", verbose, ...)
        }

    }

    for(t in seq_along(config$years)){

        if(config$years[t] > 2017){
            ppl_estimate_daoo(sp_type, year = config$years[t], config, term = "short", verbose, ...)
        }

    }
}
