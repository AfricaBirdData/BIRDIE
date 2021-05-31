#' Prepare CWAC data to fit a state-space model
#'
#' @param counts A data frame with raw CWAC count data
#' @param species An optional vector of species codes. Only records for these species are prepared for fitting an SSM, others are discarded.
#'
#' @return A tibble with clean and prepared data for fitting an SSM (e.g. filled with missing years)
#' @export
#'
#' @examples
#' counts <- getCwacSiteCounts(26352535)
#' prepSsmData(counts)
#' prepSsmData(counts, species = 212)
#' prepSsmData(counts, species = c(212, 50))
prepSsmData <- function(counts, species = NULL){

    if(!is.null(species)){
        counts <- dplyr::filter(counts, spp %in% species)
    }

    if(length(unique(counts$spp)) == 1){
        sp_code <- unique(counts$spp)
        sp_name <- paste(unique(counts$taxon.Common_species), unique(counts$taxon.Common_group))
    } else {
        sp_name <- "multi"
    }

    # Prepare season info -----------------------------------------------------

    # Some records are classified as "O" (other). We are only interested in summer and winter
    # so we filter out "O"
    counts <- counts %>%
        dplyr::filter(Season != "O")

    # We also create a numeric season variable for season
    counts <- counts %>%
        dplyr::mutate(season_id = ifelse(Season == "W", 2, 1),
               # and a year variable
               year = lubridate::year(startDate))


    # Calculate seasonal counts -----------------------------------------------

    # Sum all counts per card
    counts <- counts %>%
        dplyr::group_by(card, year, season_id) %>%
        dplyr::summarize(count = sum(count)) %>%
        dplyr::ungroup()

    # Are there more than one card per season and year?
    # It doesn't seem to, but we should be careful about this
    if(identical(length(unique(counts$year)),
                 counts %>%
                 dplyr::group_by(year, season_id) %>%
                 dplyr::summarize(n = dplyr::n()) %>%
                 dplyr::pull(n) %>%
                 sum())){
        warning("More than one card per year. Data is not appropiate for SSM")
    }

    # Fill in missing years
    counts <- counts %>%
        tidyr::complete(year = min(year):max(year), season_id)

    # Add a column with species name
    counts$spp <- sp_name

    return(counts)
}
