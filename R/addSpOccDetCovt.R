#' Add detection covariate to spOccupancy data list
#'
#' @param spOcc_data An spOccupancy data list
#' @param covt_data A data frame with the ABAP visits used to create spOcc_data.
#' The data frame must contain the columns 'Pentad', 'StartDate' and the covariate
#' that we want to add to spOcc_data. If the data is multi-season, then we also
#' need the variable that identifies the season. Only one covariate at a time
#' is allowed at the moment.
#' @param seasons The name of the variable used to identify the seasons in
#' multi-season data sets. Defaults to NULL.
#'
#' @return An spOccupancy data list with the additional detection covariate.
#' @export
#'
#' @examples
addSpOccDetCovt <- function(spOcc_data, covt_data, seasons = NULL){

    # Sort data frame for consistency with other functions
    covt_data <- covt_data %>%
        dplyr::arrange(Pentad, StartDate)

    # Extract unique pentads
    pentad_id <- dimnames(spOcc_data$y)[[1]]
    n_sites <- length(pentad_id)

    # Extract maximum number of visits in a single season
    if(is.null(seasons)){

        max_visits <- dim(spOcc_data$y)[2]
        n_seasons <- 1
        season_vec <- 1

        covt_data <- covt_data %>%
            dplyr::mutate(season = as.character(1))

    } else {

        if(!seasons %in% names(covt_data)){
            stop(paste("The field", seasons, "specified as season, is not present in the data"))
        }

        # Extract seasons and calculate max number of visits per season
        # Note that we will need to match to dimnames of spOcc_data with are characters
        covt_data <- covt_data %>%
            dplyr::mutate(season = !!as.name(seasons)) %>%
            dplyr::mutate(season = as.character(season))

        max_visits <- dim(spOcc_data$y)[3]
        n_seasons <- dim(spOcc_data$y)[2]
        season_vec <- dimnames(spOcc_data$y)[[2]]

        # Sometimes there might not be names. We make new default ones
        if(is.null(season_vec)){
            season_vec <- seq_len(n_seasons)
        }

    }

    # Aux padding vector
    vpad <- rep(NA, max_visits)

    # Aux pentad column to join yearly results
    pentad_col <- data.frame(Pentad = pentad_id)

    # Extract variable name
    var_names <- names(covt_data)[which(!names(covt_data) %in% c("Pentad", "StartDate", "season", seasons))]

    for(v in seq_along(var_names)){

        var_name <- var_names[v]

        ## Create empty array
        covt_out <-  array(NA, dim = c(n_sites, n_seasons, max_visits),
                           dimnames = list(pentad_id, season_vec, seq_len(max_visits)))

        for(k in seq_along(season_vec)){

            ## Create dataframe to format
            season_data <- covt_data %>%
                dplyr::filter(season == season_vec[k]) %>%
                # dplyr::select(Pentad, Spp, TotalHours, StartDate) %>%
                # dplyr::mutate(Spp = ifelse(Spp == "-", 0L, 1L),
                #               julian_day = lubridate::yday(StartDate)) %>%
                dplyr::right_join(pentad_col, by = "Pentad") %>%      # Join in Pentads missing for the year
                dplyr::nest_by(Pentad) %>%
                dplyr::arrange(Pentad) %>%
                dplyr::mutate(varpad = list(head(c(data[, var_name, drop = TRUE], vpad), max_visits)))

            ## Extract total hours
            covt_out[,k,] <- do.call("rbind", season_data$varpad)

        }

        # If there are no seasons then return a matrix
        if(length(season_vec) == 1){

            covt_out <- as.matrix(covt_out[, 1, ])

        }


        # Make data list
        if(is.null(spOcc_data$det.covs)) {
            spOcc_data$det.covs <- covt_out
        } else {
            spOcc_data$det.covs <- c(spOcc_data$det.covs, list(covt_out))
            names(spOcc_data$det.covs)[length(spOcc_data$det.covs)] <- var_name
        }
    }
    return(spOcc_data)

}
