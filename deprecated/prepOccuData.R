#' Prepare Spatial Occupancy Data
#'
#' @description Prepares data to fit a spatial occupancy model using the stocc
#' package.
#' @param .occudata A dataframe with at least seven columns: "Pentad", "lon",
#' "lat", "StartDate", "TotalHours", "detc" (1 for detection, 0 for
#' non-detection).
#' @param .year The year that we want to analyse
#' @param .covts A data frame with as many rows as sites (Pentads) and
#' site-level covariates that we want to use for fitting the model. The
#' covariates must have the year they correspond to in their name.
#'
#' @return An object of class so.data that we can feed into the function
#' stocc::spatial.occupancy.
#' @export
#'
#' @examples
prepOccuData <- function(.occudata, .year, .covts){

    # Remove geometry from covariates if it exists
    if("sf" %in% class(.covts)){
        .covts <- sf::st_drop_geometry(.covts)
    }

    # Prepare visit data
    visit_data <- .occudata %>%
        dplyr::filter(year %in% .year) %>%
        dplyr::select(Pentad, lon, lat, StartDate, TotalHours, detc) %>%
        dplyr::mutate(year = lubridate::year(StartDate),
                      month = lubridate::month(StartDate)) %>%
        as.data.frame()

    # Add basis values of non-linear effect for time of year
    spl_bs <- mgcv::cSplineDes(x = visit_data$month,
                               knots = seq(1, 12, length.out = 12),
                               ord = 4, derivs = 0)
    colnames(spl_bs) <- paste0("cyclic.", 1:ncol(spl_bs))
    visit_data <- cbind(visit_data, spl_bs)

    # Prepare site data
    site_data <- .covts %>%
        dplyr::select(Pentad = Name, lon, lat,
                      dplyr::contains(as.character(.year))) %>%
        dplyr::mutate(intcp = 1) %>%
        as.data.frame()

    so_data <- stocc::make.so.data(visit.data = visit_data,
                                   site.data = site_data,
                                   names = list(visit = list(site = "Pentad", obs = "detc"),
                                                site = list(site = "Pentad", coords = c("lon", "lat"))))

    return(so_data)

}
