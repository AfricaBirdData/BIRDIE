#' Annotate site data with Google Earth Engine with a time limit
#'
#' @param ee_collection A character string with the name of a GEE collection.
#' @param band The band of the collection we are going to use to annotate our data
#' @param last_year The last year provided by the GEE image collection
#' @param reducer The reducer we will use to summarize the images of the image collection
#' @param unmask GEE masks missing values, which means they are not used for
#' computing means, counts, etc. Sometimes we might want to avoid this behaviour
#' and use 0 instead of NA. If so, set unmask to TRUE.
#' @param monitor Logical. If TRUE (default) monitoring messages produced
#' by `rgee` will displayed. If FALSE, only high-level messages will be displayed.
#' @param config A list with pipeline configuration parameters.
#' See \code{\link{configPipeline}}
#' @details
#' The function annotates data with the corresponding year if it is present in the
#' GEE image collection. Years after the last year present in the config object are annotated
#' with the last year present in the image collection.
#' @return
#' @export
#'
#' @examples
addVarEETimeLimit <- function(ee_collection, band, last_year, reducer, unmask,
                              monitor, config){

    message(paste("Using", last_year, "for", band, "for years after", last_year))

    years_after <- config$years[config$years > last_year]
    years_before <- config$years[!config$years > last_year]

    # Add year "last year" to retrieve info from the last year available
    if(length(years_before) == 0){
        years_before <- last_year
    }

    if(length(years_before) == 1){

        dates <- paste0(range(years_before) + c(0,1), "-01-01")

        ee_layer <- ee_collection$
            select(band)$
            filterDate(dates[1], dates[2])$
            first()


        # Annotate with image
        out <- ABDtools::addVarEEimage(ee_feats = ee_pentads,
                                       image = ee_layer,
                                       reducer = reducer,
                                       bands = band,
                                       monitor = monitor,
                                       unmask = unmask)

        # Rename because variables are named after reducer when annotated with
        # images and not with year, like when annotated with collections
        var_name <- paste0(band, "_", years_before)

        out <- out %>%
            dplyr::rename(!!var_name := get(paste0(band, "_", reducer)))

    } else {

        yr_range <- c(min(years_before), max(years_before))

        # Make multiband image
        stackCollection <- ABDtools::EEcollectionToMultiband(collection = ee_collection,
                                                             dates = paste0(yr_range + c(0,1), "-01-01"),
                                                             band = band,
                                                             group_type = "year",
                                                             groups = years_before,
                                                             reducer = reducer,
                                                             unmask = unmask)

        # Annotate with image
        out <- ABDtools::addVarEEimage(ee_pentads, stackCollection, reducer, monitor = monitor)

    }

    # Create columns for years > last_year
    for(y in seq_along(years_after)){

        var_name <- paste0(band, "_", years_after[y])

        out <- out %>%
            dplyr::mutate(!!var_name := get(paste0(band, "_", last_year)))

    }

    return(out)

}
