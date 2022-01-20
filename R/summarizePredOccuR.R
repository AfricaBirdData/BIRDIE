#' Summarize occuR predictions for sites
#'
#' @description occuR produces estimates of detection probability for each
#' visit and occupancy probabilities for each site. This functions takes these
#' predictions and creates detection probabilities, occupancy probabilities and
#' realized occupancy probabilities (occupancy probabilities conditional on
#' observed data) for each site and occasion.
#' @param pred_p Detectopm probabilities estimated from an occuR model. If it is
#' a matrix it is interpreted as bootstrap predictions, with one row for each
#' bootstrap sample. The argument 'quants' must then be provided.
#' @param pred_psi Occupancy probabilities estimated from an occuR model. If it is
#' a matrix it is interpreted as bootstrap predictions, with one row for each
#' bootstrap sample. The argument 'quants' must then be provided.
#' @param pred_data A data.frame/data.table with the site data used to generate
#' occuR model predictions.
#' @param visit_data A data.frame/data.table with the visit data used to generate
#' occuR model predictions.
#' @param quants Quantiles to summarize predictions when bootstrap predictions
#' are supplied in pred_p and pred_psi. Passed as c("lower", "med", "upper").
#'
#' @return A tibble with estimates and/or quantiles for each pentad in site_data:
#' - psi: occupancy probability,
#' - p: detection probability,
#' - occu: realized occupancy (occupancy conditional on data).
#' @export
#'
#' @examples
summarizePredOccuR <- function(pred_p, pred_psi, pred_data, visit_data,
                               quants = NULL){

    if(is.vector(pred_p) & is.vector(pred_psi)){

        # Calculate probability of non-detections for each pentad visited
        p_nondet <- visit_data %>%
            dplyr::mutate(p = as.vector(pred_p)) %>%
            dplyr::group_by(site, occasion) %>%
            dplyr::summarize(q = prod(1-p),
                             obs = max(obs))

        # From probability of non-detection calculate the conditional occupancy
        # probability
        pred_occu <- pred_data %>%
            as.data.frame() %>%
            dplyr::mutate(psi = as.vector(pred_psi)) %>%
            dplyr::left_join(p_nondet, by = c("site", "occasion")) %>%
            dplyr::mutate(real_occu = dplyr::case_when(obs == 1 ~ 1,
                                                       is.na(obs) ~ psi,
                                                       obs == 0 ~ psi*q / (1 - psi + psi*q)),
                          p = 1 - q) %>%
            dplyr::select(Name, year, site, occasion, psi, p, real_occu)

    } else if(is.matrix(pred_p) & is.matrix(pred_psi)){

        if(nrow(pred_p) != nrow(pred_psi)){
            stop("Different number of bootstrap samples for detection and occupancy probabilities are not allowed")
        }

        if(is.null(quants)){
            stop("There are multiple predictions for p and/or psi, so quants cannot be NULL")
        }

        # Estimate realized occupancy ---------------------------------------------

        # Calculate probability of non-detections for each pentad visited
        p_nondet <- visit_data %>%
            dplyr::mutate(ub = apply(pred_p, 2, quantile, quants[3]),
                          lb = apply(pred_p, 2, quantile, quants[1]),
                          med = apply(pred_p, 2, quantile, quants[2]),
                          est = apply(pred_p, 2, mean)) %>%
            tidyr::pivot_longer(cols = c("ub", "lb", "med", "est"),
                                names_to = "lim", values_to = "p") %>%
            dplyr::group_by(site, occasion, lim) %>%
            dplyr::summarize(q = prod(1-p),
                             obs = max(obs))

        # From probability of non-detection calculate the conditional occupancy probs
        # and plot
        pred_occu <- pred_data %>%
            as.data.frame() %>%
            dplyr::mutate(ub = apply(pred_psi, 2, quantile, quants[3]),
                          lb = apply(pred_psi, 2, quantile, quants[1]),
                          med = apply(pred_psi, 2, quantile, quants[2]),
                          est = apply(pred_psi, 2, mean)) %>%
            tidyr::pivot_longer(cols = c("ub", "lb", "med", "est"),
                                names_to = "lim", values_to = "psi") %>%
            dplyr::left_join(p_nondet, by = c("site", "occasion", "lim")) %>%
            dplyr::mutate(real_occu = dplyr::case_when(obs == 1 ~ 1,
                                                       is.na(obs) ~ psi,
                                                       obs == 0 ~ psi*q / (1 - psi + psi*q)),
                          p = 1 - q) %>%
            dplyr::select(Name, year, site, occasion, psi, p, real_occu, lim)
    } else {
        stop("pred_p and pred_psi must be either both vectors or both matrices")
    }

    return(pred_occu)

}
