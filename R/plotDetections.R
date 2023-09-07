#' Plot detection from occupancy data
#'
#' @param site_data An sf object with the sites visited in a given period.
#' @param visit_data A data frame with the visits that have occurred in that a
#' given period. visit_data must have the same name as site_data for the sites
#' and the variable must be called 'site'. Detections must be represented by a
#' variable named 'obs', which should be 1 when the species was detected and 0
#' otherwise.
#'
#' @return
#' @export
#'
#' @examples
plotDetections <- function(site_data, visit_data){

    print(
        site_data %>%
            dplyr::left_join(visit_data %>%
                                 dplyr::group_by(site) %>%
                                 dplyr::summarise(everobs = sum(obs)) %>%
                                 dplyr::mutate(everobs = ifelse(everobs > 0, 1, 0)),
                             by = 'site') %>%
            ggplot2::ggplot() +
            ggplot2::geom_sf(ggplot2::aes(fill = factor(everobs)), lwd = 0.01) +
            ggplot2::scale_fill_viridis_d(name = "Obs", na.value = "grey")
    )

}
