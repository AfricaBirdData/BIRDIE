#' Predict from spOccupancy model fit
#'
#' @inheritParams ppl_summarise_occu
#'
#' @return A list with two elements: 1) posterior occupancy probability samples for
#' the South African pentads, and 2) posterior detection probability samples for each
#' visit in the ABAP data for the year `year_sel`.
#' @export
#'
#' @examples
predictSpOccu <- function(fit, sp_code, year_sel, config, ...){

    varargs <- list(...)

    # Prepare prediction data -------------------------------------------------

    # Load site data for prediction
    sitefile <- file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv"))
    visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

    # Load data
    sitedata <- utils::read.csv(sitefile, check.names = FALSE)
    visitdata <- utils::read.csv(visitfile, check.names = FALSE)
    detdata <- utils::read.csv(detfile, check.names = FALSE)


    # Format for occupancy modelling ------------------------------------------

    occudata <- BIRDIE::createOccuData(sp_code = sp_code,
                                       years = year_sel,
                                       site_data = sitedata,
                                       visit_data = NULL,
                                       config = config,
                                       force_abap_dwld = FALSE)

    # Add detection info to visit data and subset years
    occudata$visit <- visitdata %>%
        dplyr::left_join(detdata,
                         by = c("CardNo", "StartDate", "year", "Pentad")) %>%
        dplyr::filter(year == year_sel)

    # Save codes of pentads we are predicting occupancy and detection for
    occ_pentads <- occudata$site$Pentad
    det_pentads <- occudata$visit$Pentad

    # Save also detections in visits
    obs_visit <- occudata$visit$obs

    # Select occupancy variables to be included in the model
    tt_occ <- stats::terms(stats::reformulate(config$occ_mod))
    occ_vars <- attr(tt_occ, "term.labels")
    occ_vars <- gsub(".* \\| ", "", occ_vars)

    occudata$site <- occudata$site %>%
        dplyr::select(dplyr::all_of(occ_vars))

    # Select detection variables to be included in the model
    tt_det <- stats::terms(stats::reformulate(config$det_mod))
    det_vars <- attr(tt_det, "term.labels")
    det_vars <- gsub(".* \\| ", "", det_vars)

    occudata$visit <- occudata$visit %>%
        dplyr::select(dplyr::all_of(det_vars))

    # Scale and center as for model
    scale_fct_occ <- fit$occ.scale
    occudata$site[,names(scale_fct_occ$center)] <- scale(occudata$site[,names(scale_fct_occ$center)],
                                                         center = scale_fct_occ$center,
                                                         scale =  scale_fct_occ$scale)

    scale_fct_det <- fit$det.scale
    occudata$visit[,names(scale_fct_det$center)] <- scale(occudata$visit[,names(scale_fct_det$center)],
                                                          center = scale_fct_det$center,
                                                          scale =  scale_fct_det$scale)

    # Add intercept
    occudata <- purrr::map(occudata, ~ .x %>%
                               dplyr::mutate(intcp = 1) %>%
                               dplyr::select(intcp, dplyr::everything()))


    pred_data <- list(psi = NA, p = NA)
    pred_data$psi <- spOccupancy:::predict.PGOcc(fit, as.matrix(occudata$site), ignore.RE = FALSE, type = "occupancy")$psi.0.samples
    pred_data$p <- spOccupancy:::predict.PGOcc(fit, as.matrix(occudata$visit), ignore.RE = FALSE, type = "detection")$p.0.samples

    # Add pentad information
    attr(pred_data$psi, "pentads") <- occ_pentads
    attr(pred_data$p, "pentads") <- det_pentads

    # And detections information
    attr(pred_data$p, "obs") <- obs_visit

    # Add year information
    attr(pred_data$psi, "year") <- year_sel
    attr(pred_data$p, "year") <- year_sel

    return(pred_data)

}
