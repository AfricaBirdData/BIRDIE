#' Prepare CWAC data to fit a state-space model
#'
#' @param counts A data frame with raw CWAC count data
#' @param spp_sel An optional vector of species codes. Only records for these species are prepared for fitting an SSM, others are discarded.
#'
#' @return A tibble with clean and prepared data for fitting an SSM (e.g. filled with missing years)
#' @export
#'
#' @examples
#' counts <- barberspan
#' prepSsmData(counts)
#' prepSsmData(counts, spp_sel = 212)
#' prepSsmData(counts, spp_sel = c(212, 50))
prepSsmData <- function(counts, spp_sel = NULL){

    # Prepare variable names
    counts <- counts %>%
        dplyr::rename(spp = SppRef) %>%
        dplyr::rename_with(.fn = ~tolower(.x))

    # Prepare species output name
    if(length(spp_sel) == 1){
        sp_name <- paste(unique(counts[counts$spp == spp_sel, "common_species"]),
                         unique(counts[counts$spp == spp_sel, "common_group"]))
    } else {
        sp_name <- "multi"
    }

    # Prepare season info -----------------------------------------------------

    # Some records are classified as "O" (other). We are only interested in summer and winter
    # so we filter out "O"
    counts <- counts %>%
        dplyr::filter(season != "O")

    # We also create a numeric season variable for season
    counts <- counts %>%
        dplyr::mutate(season_id = ifelse(season == "W", 2, 1),
                      # and a year variable
                      year = lubridate::year(startdate))

    # Are there more than one card per season and year?
    # It doesn't seem to, but we should be careful about this
    if(1 < counts %>%
       dplyr::group_by(year, season_id) %>%
       dplyr::summarize(n = dplyr::n_distinct(card)) %>%
       dplyr::pull(n) %>%
       max()){
        warning("More than one card per year. Data is not appropiate for SSM")
    }

    # Calculate seasonal counts

    # Sum all counts per card
    counts_all_spp <- counts %>%
        dplyr::group_by(card, year, season_id) %>%
        dplyr::summarize(count = sum(count)) %>%
        dplyr::ungroup() %>%
        tidyr::complete(year = min(year):max(year), season_id) # Fill in missing years


    # Filter species ----------------------------------------------------------

    if(!is.null(spp_sel)){
        counts_sp <- counts %>%
            dplyr::filter(spp %in% spp_sel) %>%
            dplyr::group_by(card, year, season_id) %>%
            dplyr::summarize(count = sum(count)) %>%
            dplyr::ungroup()

        counts <- counts_all_spp %>%
            dplyr::left_join(counts_sp[,c("card", "count")],
                             by = c("card")) %>%
            dplyr::mutate(count.y = ifelse(is.na(count.y) & !is.na(count.x), 0, count.y)) %>%
            dplyr::rename(count = count.y) %>%
            dplyr::select(-count.x)

    } else {
        counts <- counts_all_spp
    }

    # Add a column with species name
    counts$spp <- sp_name

    return(counts)
}
