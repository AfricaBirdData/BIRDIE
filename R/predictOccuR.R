#' Predict occupancy from occuR model
#'
#' @param fit A occuR model fit. See \link[occuR]{fit_occu}
#' @param occuRdata A list produced by \link{prepDataOccuR}
#' @param pred_sites A spatial `sf` object with the sites to predict occupancy
#' for.
#' @param years A vector with years to predict for.
#' @param scaling Logical; whether it is necessary to scale the variables
#' before predicting. Scaling factors will be extracted from the attributes of
#' occuRdata. See \link{scale}.
#' @param quants A vector of three quantiles of occupancy predictions to be
#' computed. Passed as c("lower", "med", "upper").
#' @param boots Number of bootstrap samples to calculate quantiles.
#'
#' @return A tibble with estimates and quantiles for each pentad in pred_sites:
#' - psi: occupancy probability,
#' - p: detection probability,
#' - occu: realized occupancy (occupancy conditional on data).
#' @export
#'
#' @examples
predictOccuR <- function(fit, occuRdata, pred_sites, years, scaling = FALSE,
                         quants = c(0.025, 0.5, 0.975), boots = 1000){

    # Prepare data to predict psi and p
    gm <- pred_sites %>%
        dplyr::select(Name)

    # Separate variables into columns and add necessary covariates
    pred_data <- pred_sites %>%
        sf::st_drop_geometry()  %>%
        dplyr::group_by(Name) %>%
        dplyr::mutate(site = dplyr::cur_group_id()) %>%
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
        dplyr::group_by(year) %>%
        dplyr::mutate(occasion = cur_group_id()) %>%
        dplyr::ungroup() %>%
        data.table::as.data.table()

    # Predict
    pred <- predict(fit, occuRdata$visit,  pred_data, nboot = boots)


    # Estimate realized occupancy ---------------------------------------------

    # Calculate probability of non-detections for each pentad visited
    p_nondet <- occuRdata$visit %>%
        dplyr::select(Pentad, site, occasion, obs) %>%
        dplyr::mutate(ub = apply(pred$pboot, 2, quantile, quants[3]),
                      lb = apply(pred$pboot, 2, quantile, quants[1]),
                      med = apply(pred$pboot, 2, quantile, quants[2]),
                      est = pred$p) %>%
        tidyr::pivot_longer(cols = -c(Pentad, site, occasion, obs),
                            names_to = "lim", values_to = "p") %>%
        dplyr::group_by(Pentad, site, occasion, lim) %>%
        dplyr::summarize(q = prod(1-p),
                         obs = max(obs))

    # From probability of non-detection calculate the conditional occupancy probs
    # and plot
    pred_occu <- pred_data %>%
        as.data.frame() %>%
        dplyr::left_join(gm, by = "Name") %>%
        sf::st_sf() %>%
        dplyr::mutate(ub = apply(pred$psiboot, 2, quantile, quants[3]),
                      lb = apply(pred$psiboot, 2, quantile, quants[1]),
                      med = apply(pred$psiboot, 2, quantile, quants[2]),
                      est = pred$psi[,1]) %>%
        tidyr::pivot_longer(cols = c("ub", "lb", "med", "est"),
                            names_to = "lim", values_to = "psi") %>%
        dplyr::left_join(p_nondet, by = c("Name" = "Pentad", "occasion", "lim")) %>%
        dplyr::mutate(real_occu = dplyr::case_when(obs == 1 ~ 1,
                                                   is.na(obs) ~ psi,
                                                   obs == 0 ~ psi*q / (1 - psi + psi*q)),
                      p = 1 - q) %>%
        dplyr::select(Name, year, psi, p, real_occu, lim)

    return(pred_occu)
}
