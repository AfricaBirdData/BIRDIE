#' Prepare CWAC data to fit a continuous-time state-space model
#'
#' @param counts A data frame with raw CWAC count data
#' @param species An optional vector of species codes. Only records for these
#' species are prepared for fitting an SSM, others are discarded.
#'
#' @return A tibble with clean and prepared data for fitting a continuous-time
#' SSM (e.g. filled with missing years)
#' @export
#'
#' @examples
#' counts <- barberspan
#' prepCtSsmData(counts)
#' prepCtSsmData(counts, species = 212)
#' prepCtSsmData(counts, species = c(212, 50))
prepCtSsmData <- function(counts, species = NULL){

    # Prepare species output name
    if(length(species) == 1){
        sp_name <- paste(unique(counts[counts$spp == species, "taxon.Common_species"]),
                         unique(counts[counts$spp == species, "taxon.Common_group"]))
    } else {
        sp_name <- "multi"
    }

    # Prepare season info -----------------------------------------------------

    # Prepare times and time increments
    datetimes <- data.frame(startDate = unique(counts$startDate)) %>%
        dplyr::arrange(startDate) %>%
        dplyr::mutate(dt = as.numeric(difftime(dplyr::lead(startDate), startDate, units = "days")),
                      # dt = ifelse(is.na(dt), 0, dt),
                      t = c(0, cumsum(dt[!is.na(dt)])))

    counts <- counts %>%
        dplyr::left_join(datetimes, by = "startDate")

    # Sum all counts per card
    cc <- counts %>%
        dplyr::group_by(t) %>%
        dplyr::summarize(card = unique(card),
                  date = unique(startDate),
                  season = unique(Season),
                  dt = unique(dt),
                  count = sum(count)) %>%
        dplyr::mutate(year = lubridate::year(date))

    # Create dseason, number of seasons between observations
    dseason <- vector("integer", length = nrow(cc))
    ss <- c(diff(cc$year[cc$season == "S"]), NA)
    ww <- c(diff(cc$year[cc$season == "W"]), NA)
    ids <- 1
    idw <- 1

    for(i in 1:nrow(cc)){

        # update dseason
        if(cc$season[i] == "S"){
            dseason[i] <- ss[ids]
            ids <- ids + 1
        }

        if(cc$season[i] == "W"){
            dseason[i] <- ww[idw]
            idw <- idw + 1
        }

    }

    cc <- cc %>%
        dplyr::mutate(dseason = dseason)

    # We also create dummy variables for the seasons
    cc <- cc %>%
        dplyr::mutate(summer = ifelse(season == "S", 1, 0),
                      winter = ifelse(season == "W", 1, 0))


    # Filter species ----------------------------------------------------------

    if(!is.null(species)){
        counts_sp <- counts %>%
            dplyr::filter(spp %in% species) %>%
            dplyr::group_by(card) %>%
            dplyr::summarize(count = sum(count)) %>%
            dplyr::ungroup()

        counts <- cc %>%
            dplyr::left_join(counts_sp[,c("card", "count")],
                             by = c("card")) %>%
            dplyr::mutate(count.y = ifelse(is.na(count.y) & !is.na(count.x), 0, count.y)) %>%
            dplyr::rename(count = count.y) %>%
            dplyr::select(-count.x)

    } else {
        counts <- cc
    }

    # Add a column with species name
    counts$spp <- sp_name

    # Add 0.1 to zero counts otherwise log goes to minus infinity
    counts <- counts %>%
        dplyr::mutate(count = count)

    return(counts)
}
