#' Add occupancy visit covariates
#'
#' @param visit_data A dataframe with occupancy data
#' @param sites A spatial object with the sites where visit_data was collected.
#' @param covts A character vector with the names of the covariates to be
#' extracted. The names must correspond to the names contained the covariate
#' file name (e.g. "prcp", "tmax", etc.). At the moment the function assumes
#' that covariates are found at: "covts_dir_file_fix[1]covts_file_fix[2].rds"
#' @param covts_dir Directory where covariates are found.
#' @param file_fix A character vector with two elements correspondint to
#' additional characters that are found before and after covts.
#'
#' @returnThe visit_data dataframe with additional columns corresponding to
#' the new variables.
#' @export
#'
#' @examples
addOccVisitCovts <- function(visit_data, sites, covts, covts_dir, file_fix){

    # Identify all different year-month combinations
    visit_data <- visit_data %>%
        dplyr::mutate(month = dplyr::if_else(month < 10,
                                      paste0("0", month),
                                      as.character(month))) %>%
        dplyr::mutate(yearmonth = paste0("X", year, ".", month))

    yearmonth <- visit_data %>%
        dplyr::distinct(yearmonth) %>%
        dplyr::pull(yearmonth)

    visit_data <- dplyr::nest_by(visit_data, yearmonth)

    for(i in seq_along(covts)){

        # Load covariate
        covt <- readRDS(paste0(covts_dir, file_fix[1], covts[i],
                               file_fix[2], ".rds"))

        # Subset and annotate data
        print(paste("Extracting variable", covts[i]))
        p <- progressr::progressor(steps = nrow(visit_data))

        visit_data$data <- furrr::future_map2(
            yearmonth, visit_data$data,
            ~extractVisitCovt(.sub_by = .x, .visit_data = .y, .sites = sites,
                             .covt = covt, .name = covts[i], .p = p))
    }

    # Unnest
    visit_data <- visit_data %>%
        dplyr::summarize(data)

    return(visit_data)

}

