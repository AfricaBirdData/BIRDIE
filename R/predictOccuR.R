#' Predict from occuR model fit
#'
#' @inheritParams ppl_summarise_occu
#'
#' @return A list with two elements: 1) posterior occupancy probability samples for
#' the South African pentads, and 2) posterior detection probability samples for each
#' visit in the ABAP data for the year `year_sel`.
#'
#' @return
#' @export
#'
#' @examples
predictOccuR <- function(fit, sp_code, year, config, ...){

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
        dplyr::filter(year == year_sel) %>%
        dplyr::mutate(Spp = ifelse(obs == 0, "-", "1"))

    # # Save codes of pentads we are predicting occupancy and detection for
    # occ_pentads <- occudata$site$Pentad
    # det_pentads <- occudata$visit$Pentad
    #
    # # Save also detections in visits
    # obs_visit <- occudata$visit$obs

    # Format to occuR
    occuRdata <- prepOccuRData(occudata$site, occudata$visit, config, spatial = FALSE,
                               sp_sites = NULL, scale = FALSE, keep_sites = TRUE)

    # # Select occupancy variables to be included in the model
    # tt_occ <- stats::terms(stats::reformulate(config$occ_mod))
    # occ_vars <- attr(tt_occ, "term.labels")
    # occ_vars <- gsub(".* \\| ", "", occ_vars)
    #
    # occuRdata$site <- occuRdata$site %>%
    #     dplyr::select(dplyr::all_of(occ_vars))
    #
    # # Select detection variables to be included in the model
    # tt_det <- stats::terms(stats::reformulate(config$det_mod))
    # det_vars <- attr(tt_det, "term.labels")
    # det_vars <- gsub(".* \\| ", "", det_vars)
    #
    # occuRdata$visit <- occuRdata$visit %>%
    #     dplyr::select(dplyr::all_of(det_vars))

    # Scale and center as for model
    occuRdata <- lapply(occuRdata, as.data.frame)
    scale_fct_occ <- fit$occ.scale
    scale_fct_occ <- lapply(scale_fct_occ, function(x) x[-grep(":", names(x))])   # For occuR we need to remove interactions
    occuRdata$site[,names(scale_fct_occ$center)] <- scale(occuRdata$site[,names(scale_fct_occ$center)],
                                                         center = scale_fct_occ$center,
                                                         scale =  scale_fct_occ$scale)

    scale_fct_det <- fit$det.scale
    scale_fct_det <- lapply(scale_fct_det, function(x) x[-grep(":", names(x))])   # For occuR we need to remove interactions
    occuRdata$visit[,names(scale_fct_det$center)] <- scale(occuRdata$visit[,names(scale_fct_det$center)],
                                                          center = scale_fct_det$center,
                                                          scale =  scale_fct_det$scale)

    occuRdata <- lapply(occuRdata, data.table::as.data.table)

    occuRdata$site <- occuRdata$site %>%
        dplyr::mutate(site_id = factor(site_id))
    occuRdata$visit <- occuRdata$visit %>%
        dplyr::mutate(obs_id = factor(obs_id),
                      site_id = factor(site_id))


    tryCatch({

        # Predict
        pred_occu <- predict(fit, occuRdata$visit, occuRdata$site,
                             include_re = TRUE, new_levels = TRUE, nboot = 1000)

    }, error = function(e){
        sink(file.path(config$out_dir, sp_code, paste0("failed_pred_", config$package, "_", sp_code,"_", year, ".txt")))
        print(e)
        sink()}) # TryCatch predict

    # Add pentad information
    # attr(pred_data$psi, "pentads") <- occ_pentads
    # attr(pred_data$p, "pentads") <- det_pentads

    # And detections information
    # attr(pred_data$p, "obs") <- obs_visit

    # Add year information
    # attr(pred_data$psi, "year") <- year_sel
    # attr(pred_data$p, "year") <- year_sel

    attr(pred_occu$psi, "pentad") <- occuRdata$site$pentad
    attr(pred_occu$psi, "year") <- occuRdata$site$year
    attr(pred_occu$p, "pentad") <- occuRdata$visit$pentad
    attr(pred_occu$p, "obs") <- occuRdata$visit$obs
    attr(pred_occu$psiboot, "pentad") <- occuRdata$site$pentad
    attr(pred_occu$psiboot, "year") <- occuRdata$site$year
    attr(pred_occu$pboot, "pentad") <- occuRdata$visit$pentad
    attr(pred_occu$pboot, "obs") <- occuRdata$visit$obs

    return(pred = pred_occu)

}
