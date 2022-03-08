#' Prepare CWAC data to fit a state-space model
#'
#' @param counts A data frame with raw CWAC count data
#' @param spp_sel An optional vector of species codes. Only records for these
#' species are prepared for fitting an SSM, others are discarded.
#' @param keep A vector of variables to keep after the processing. Useful, if
#' there are covariates of interest. If NULL, only year, season, start date,
#' count and species name are returned.
#'
#' @return A tibble with clean and prepared data for fitting an SSM (e.g. filled with missing years)
#' @export
#'
#' @examples
#' counts <- barberspan
#' prepSsmData(counts)
#' prepSsmData(counts, spp_sel = 212)
#' prepSsmData(counts, spp_sel = c(212, 50))
prepSsmData <- function(counts, spp_sel = NULL, keep = NULL){

    # Save dataframe to save all the original variables
    counts_orig <- counts

    # Prepare variables and sort data by date
    counts <- counts %>%
        dplyr::mutate(year = lubridate::year(StartDate)) %>%
        dplyr::arrange(StartDate)

    # Prepare species output name
    if(!is.null(spp_sel) && length(spp_sel) == 1){
        sp_name <- counts %>%
            dplyr::filter(!is.na(Common_species),
                          SppRef == spp_sel) %>%
            dplyr::distinct(sp_name = paste(Common_species, Common_group)) %>%
            dplyr::pull(sp_name)
    } else {
        sp_name <- "multi"
    }


    # Check one card per day -------------------------------------------------

    if((counts %>%
        dplyr::filter(!is.na(Card)) %>%
        dplyr::count(Card, StartDate) %>%
        dplyr::count(Card) %>%
        dplyr::pull(n) %>% max() > 1)){
        stop("Some cards are repeated in different days!")
    } else if((counts %>%
               dplyr::filter(!is.na(Card)) %>%
               dplyr::count(Card, StartDate) %>%
               dplyr::count(StartDate) %>%
               dplyr::pull(n) %>% max() > 1)){
        stop("Some dates have more than one card!")
    }

    # Filter species ----------------------------------------------------------

    if(!is.null(spp_sel)){

        counts <- counts %>%
            dplyr::mutate(count = dplyr::case_when(SppRef %in% spp_sel ~ Count,
                                                   is.na(SppRef) ~ NA_integer_,
                                                   TRUE ~ 0L))

        counts_sum <- counts %>%
            dplyr::group_by(year, Season, StartDate) %>%
            dplyr::summarize(count = sum(count, na.rm = F)) %>%
            dplyr::ungroup()

    } else {

        counts_sum <- counts %>%
            dplyr::group_by(year, Season, StartDate) %>%
            dplyr::summarize(count = sum(Count, na.rm = F)) %>%
            dplyr::ungroup()

    }

    # Add a column with species name
    counts_sum$spp <- sp_name

    # Put back original variables
    if(!is.null(keep)){
        counts_sum <- counts_sum %>%
            dplyr::left_join(counts_orig %>%
                                 dplyr::select(Season, StartDate, dplyr::all_of(keep)) %>%
                                 dplyr::distinct(),
                             by = c("Season", "StartDate"))
    }

    return(counts_sum)
}
