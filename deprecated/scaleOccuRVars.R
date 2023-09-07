#' Scale variables in occuR-like data
#'
#' @param site_data A dataframe with occupancy site data.
#' @param visit_data A dataframe with occupancy visit data.
#' @param scale_vars A named list with names 'site' and 'visit'. Each element
#' must contain the variables to be scaled in the site and visit data
#' respectively. If one element is NULL, no variables will be scaled in the data
#' corresponding to that element.
#'
#' @return A named list containing scaled and centred 'site' and 'visit' data.
#' @export
#'
#' @examples
scaleOccuRVars <- function(site_data, visit_data, scale_vars){

    if(!is.null(scale_vars$visit)){
        visit_data <- visit_data %>%
            dplyr::mutate(dplyr::across(.col = dplyr::all_of(scale_vars$visit), .fns = ~scale(.x)))
        attr(visit_data, "scaling") <- TRUE
    } else {
        attr(visit_data, "scaling") <- FALSE
    }

    if(!is.null(scale_vars$site)){
        site_data <- site_data %>%
            dplyr::mutate(dplyr::across(.col = dplyr::all_of(scale_vars$site), .fns = ~scale(.x)))
        attr(site_data, "scaling") <- TRUE
    } else {
        attr(site_data, "scaling") <- FALSE
    }

    return(list(site = site_data,
                visit = visit_data))

}
