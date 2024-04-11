#' Gather years from variables with year
#'
#' @param x A dataframe with columns that represent a variable on a given year.
#' Note that there must be a unique key that identifies each case uniquely.
#' @param vars Variables to separate in variable and year (in different columns).
#' @param sep A character that identifies the separator between the name of the
#' variable and the year.
#'
#' @return A dataframe with separate columns for variables and years. There will
#' be a column that represents the year and several columns that represent the
#' variables measured across years.
#' @export
#'
#' @examples
#' df <- data.frame(id = 1:10,
#'                  v1_2010 = rnorm(10, 0, 1),
#'                  v1_2011 = rnorm(10, 100, 10),
#'                  v2_2010 = runif(10, 0, 1),
#'                  v2_2011 = runif(10, 100, 110))
#'
#' gatherYearFromVars(x = df, vars = names(df)[-1], sep = "_")

gatherYearFromVars <- function(x, vars, sep){

    if("sf" %in% class(x)){
        stop("gatherYearFromVars does not take sf objects, please try as.data.frame(x)")
    }

    x <- x %>%
        tidyr::pivot_longer(cols = dplyr::all_of(vars)) %>%
        tidyr::separate(name, into = c("covt", "year"), sep = sep) %>%
        tidyr::pivot_wider(names_from = covt, values_from = value) %>%
        dplyr::mutate(year = as.integer(year))

    return(x)

}
