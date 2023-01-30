#' Summarise occuR occupancy estimates
#'
#' @description occuR produces estimates of detection probability for each
#' visit and occupancy probabilities for each site. This functions takes these
#' predictions and creates detection probabilities, occupancy probabilities and
#' realized occupancy probabilities (occupancy probabilities conditional on
#' observed data) for each site.
#' @param pred_psi Occupancy probabilities estimated from an a occuR model.
#' It must be a matrix with each row corresponding to an bootstrap sample and each column
#' corresponding to a pentad. This object must contain a "pentad" attribute indicating
#' which pentad each column correspond to, and an attribute "year" indicating what
#' year the probabilities correspond to. Outputs from \code{\link{predictOccuR}},
#' should be readily appropriate.
#' @param pred_p Detection probabilities estimated from a occuR model.
#' It must be a matrix with each row corresponding to an bootstrap sample and each column
#' corresponding to a visit. This object must contain a "pentad" attribute indicating
#' which pentad each column correspond to, an attribute "year" indicating what
#' year the probabilities correspond to, and an attribute "obs" indicating whether
#' the species was detected in the visit or not. Outputs from \code{\link{predictOccuR}},
#' should be readily appropriate.
#' @param quants Quantiles to summarise predictions distribution passed as c("lower", "med", "upper").
#'
#' @return A tibble with estimates and/or quantiles for each pentad in site_data:
#' - psi: occupancy probability,
#' - p: detection probability,
#' - occu: realized occupancy (occupancy conditional on data).
#' @export
#'
#' @examples
summariseOccuR <- function(pred_p, pred_psi, quants){

    # Estimate realized occupancy ---------------------------------------------

    # Calculate probability of non-detections for each pentad visited
    p_nondet <- data.frame(pentad = attr(pred_p, "pentad"),
                           obs = attr(pred_p, "obs")) %>%
        dplyr::mutate(ub = apply(pred_p, 2, quantile, quants[3]),
                      lb = apply(pred_p, 2, quantile, quants[1]),
                      med = apply(pred_p, 2, quantile, quants[2]),
                      est = apply(pred_p, 2, mean)) %>%
        tidyr::pivot_longer(cols = c("ub", "lb", "med", "est"),
                            names_to = "lim", values_to = "p") %>%
        dplyr::group_by(pentad, lim) %>%
        dplyr::summarise(q = prod(1-p),
                         obs = max(obs)) %>%
        dplyr::ungroup()

    # From probability of non-detection calculate the conditional occupancy probs
    # and plot
    pred_occu <- data.frame(pentad = attr(pred_psi, "pentad"),
                            year = attr(pred_psi, "year")) %>%
        dplyr::mutate(ub = apply(pred_psi, 2, quantile, quants[3]),
                      lb = apply(pred_psi, 2, quantile, quants[1]),
                      med = apply(pred_psi, 2, quantile, quants[2]),
                      est = apply(pred_psi, 2, mean)) %>%
        tidyr::pivot_longer(cols = c("ub", "lb", "med", "est"),
                            names_to = "lim", values_to = "psi") %>%
        dplyr::left_join(p_nondet %>%
                             dplyr::select(pentad, lim, q, obs), by = c("pentad", "lim")) %>%
        dplyr::mutate(real_occu = dplyr::case_when(obs == 1 ~ 1,
                                                   is.na(obs) ~ psi,
                                                   obs == 0 ~ psi*q / (1 - psi + psi*q)),
                      p = 1 - q) %>%
        dplyr::select(pentad, year, psi, p, real_occu, lim)

    return(pred_occu)

}
