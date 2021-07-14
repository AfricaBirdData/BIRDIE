#' Extract raster covariates for occupancy visit data
#'
#' @param .sub_by A character string corresponding to the name used to subset
#' the covariate (e.g. "February", "2015", etc.)
#' @param .visit_data A dataframe with occupancy visit data.
#' @param .sites A spatial object with the sites sampled to obtain .visit_data.
#' @param .covt A raster layer with the covariate to extract. If .covt is a
#' rasterStack, the layer with name equal to .sub_by with be extracted.
#' @param .name Name that will be given to the new extracted covariate
#' @param .p This is needed for progress reports. Just leave the default.
#'
#' @return The .visit_data dataframe with an additional column with the new
#' variable
#' @export
#'
#' @examples
extractVisitCovt <- function(.sub_by, .visit_data, .sites, .covt, .name, .p = p){

    # Progress report
    .p()

    sites_sel <- unique(.visit_data$SiteName)

    covt_sel <- raster::subset(.covt, .sub_by)

    sites_spt_sel <- .sites %>%
        dplyr::filter(Name %in% sites_sel)

    covt_site <- raster::extract(covt_sel, sf::as_Spatial(sites_spt_sel), fun = mean)

    sites_spt_sel <- sites_spt_sel %>%
        sf::st_drop_geometry() %>%
        dplyr::mutate(!!.name := covt_site[,1]) %>%
        dplyr::select(SiteName = Name, .name)

    .visit_data <- .visit_data %>%
        dplyr::left_join(sites_spt_sel, by = "SiteName")

    return(.visit_data)

}
