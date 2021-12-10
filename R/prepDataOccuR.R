#' Prepare data for fitting a model with occuR
#'
#' @param site_data A dataframe with occupancy site data.
#' @param visit_data A dataframe with occupancy visit data.
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

    # Keep those sites that appear in visits data and drop geometry
    site_data <- site_data %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(Name %in% sites) %>%
        dplyr::group_by(Name) %>%
        dplyr::mutate(site = dplyr::cur_group_id()) %>%
        dplyr::ungroup()

    # Separate covariates and years
    site_data <- site_data %>%
        tidyr::pivot_longer(cols = -c(Name, lon, lat, site, watocc_ever, dist_coast)) %>%
        tidyr::separate(name, into = c("covt", "year"), sep = "_") %>%
        tidyr::pivot_wider(names_from = covt, values_from = value) %>%
        dplyr::mutate(year = as.numeric(year)) %>%
        dplyr::group_by(year) %>%
        dplyr::mutate(occasion = cur_group_id()) %>%
        dplyr::ungroup()

    if(any(is.na(site_data$year))){
        warning("There might be covariates that don't change over time other than watocc_ever and dist_coast.")
    }

    # Transfer site and occasion indicators from site data
    visit_data <- visit_data %>%
        dplyr::left_join(site_data %>%
                             dplyr::select(Name, site) %>%
                             dplyr::distinct(),
                         by = c("Pentad" = "Name")) %>%
        dplyr::left_join(site_data %>%
                             dplyr::select(year, occasion) %>%
                             dplyr::distinct(),
                         by = "year")

    # Add visit indicator
    visit_data <- visit_data %>%
        dplyr::group_by(site, occasion) %>%
        dplyr::mutate(visit = row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(CardNo, Date))

    if(any(is.na(visit_data$site))){
        warning("Sites differ between site data and visit data!")
    }

    return(list(site = data.table::as.data.table(site_data),
                visit = data.table::as.data.table(visit_data)))

}
