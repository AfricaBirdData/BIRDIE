#' Prepare site data to predict from occuR occupancy model
#'
#' @param occuRdata A list produced by \link{prepDataOccuR}
#' @param pred_sites A dataframe with the sites to predict occupancy for.
#' @param years A vector with years to predict for.
#' @param scaling Logical; whether it is necessary to scale the variables
#' before predicting. Scaling factors will be extracted from the attributes of
#' occuRdata. See \link{scale}.
#'
#' @return A dataframe with site information to predict from an occuR model.
#' Sites and occasions are formatted according to occuRdata.
#' @export
#'
#' @examples
prepPredictDataOccuR <- function(occuRdata, pred_sites, years, scaling = FALSE){

    # Add site info if not present
    pred_sites <- pred_sites %>%
        dplyr::left_join(occuRdata$site %>%
                             dplyr::select(Name, site) %>%
                             dplyr::distinct(),
                         by = "Name")

    # Separate variables into columns and add necessary covariates
    pred_data <- pred_sites %>%
        tidyr::pivot_longer(cols = -c(Name, lon, lat, site, watocc_ever, dist_coast)) %>% # Note that these are hard-coded
        tidyr::separate(name, into = c("covt", "year"), sep = "_") %>%
        tidyr::pivot_wider(names_from = covt, values_from = value) %>%
        dplyr::mutate(year = as.integer(year),
                      tdiff = tmmx - tmmn) %>%
        dplyr::filter(year %in% years)

    # Scale variables
    if(scaling){
        # Extract scaling factors
        sc <- lapply(occuRdata$site, attributes)
        sc <- sc[!sapply(sc, is.null)]

        for(i in seq_along(sc)){
            pred_data[, names(sc)[i]] <- scale(x = pred_data[, names(sc)[i]],
                                               center = sc[[i]]$`scaled:center`,
                                               scale = sc[[i]]$`scaled:scale`)
        }
    }

    # Define occasion
    pred_data <- pred_data %>%
        dplyr::left_join(occuRdata$visit %>%
                             dplyr::select(year, occasion) %>%
                             dplyr::distinct(),
                         by = "year") %>%
        data.table::as.data.table()

    return(pred_data)
}
