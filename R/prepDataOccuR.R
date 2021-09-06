#' Prepare data for fitting a model with occuR
#'
#' @param site_data A dataframe with occupancy site data produced by
#' \link{prepOccSiteData}.
#' @param visit_data A dataframe with occupancy visit data produced by
#' \link{prepOccVisitData}.
#'
#' @return A list containing two data frames: one with site data and one with
#' visit data ready to use with the occuR package
#' @export
#'
#' @examples
prepDataOccuR <- function(site_data, visit_data){

    # Visited sites
    visits <- visit_data %>%
        dplyr::distinct(Pentad, year) %>%
        dplyr::mutate(keep = 1)

    sites <- unique(visit_data$Pentad)

    sitedata <- site_data %>%
        sf::st_drop_geometry()  %>%
        tidyr::pivot_longer(cols = -c(Name, lon, lat, id, water)) %>%
        tidyr::separate(name, into = c("covt", "year"), sep = "_") %>%
        tidyr::pivot_wider(names_from = covt, values_from = value) %>%
        dplyr::mutate(year = as.numeric(year)) %>%
        dplyr::left_join(visits, by = c("Name" = "Pentad", "year" = "year")) %>%
        dplyr::filter(keep == 1) %>%
        dplyr::select(-keep) %>%
        dplyr::rename(site = id) %>%
        dplyr::group_by(year) %>%
        dplyr::mutate(occasion = cur_group_id()) %>%
        dplyr::ungroup() %>%
        data.table::as.data.table()

    visitdata <- visit_data %>%
        dplyr::filter(year %in% unique(sitedata$year)) %>%
        dplyr::left_join(sitedata %>%
                             dplyr::select(Name, site) %>%
                             distinct(),
                         by = c("Pentad" = "Name")) %>%
        dplyr::left_join(sitedata %>%
                             dplyr::select(year, occasion) %>%
                             distinct(),
                         by = "year") %>%
        dplyr::group_by(site, occasion) %>%
        dplyr::mutate(visit = row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(id, CardNo, StartDate, EndDate, month_p)) %>%
        data.table::as.data.table()

    return(list(site = sitedata, visit = visitdata))

}
