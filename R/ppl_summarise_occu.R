#' Summarise predictions from occupancy model
#'
#' @inheritParams ppl_run_pipe_dst1
#' @param fit An occupancy model fit to summarise occupancy and detection
#' probabilities from.
#' @param year_sel Year in data to run model for.
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
ppl_summarise_occu <- function(fit, sp_code, sp_name, year_sel, config, ...){

    # Predict from model
    message(paste("Predicting from model", file.path(config$out_dir, sp_code, paste0("occu_fit_", year_sel, "_", sp_code, ".rds"))))
    pred_occu <- predictSpOccu(fit, sp_code, year_sel, config)

    # summarise predictions
    message("Summarising predictions")
    summ_occu <- summariseSpOcc(pred_psi = pred_occu$psi,
                                pred_p = pred_occu$p,
                                quants = c(0.025, 0.5, 0.975))

    # Save predictions

    # Subset predictions and add species
    summ_occu <- summ_occu %>%
        dplyr::mutate(species = sp_code) %>%
        dplyr::select(species, everything())

    # Save predictions
    summ_occu %>%
        write.csv(file.path(config$out_dir, sp_code, paste0("occur_pred_", year_sel, "_", sp_code, ".csv")),
                  row.names = FALSE)

    ## PLOTS

    # Download geometry if not present on disk
    pentads_file <- file.path(tempdir(), "sa_pentads.rds")

    if(file.exists(pentads_file)){
        pentads_sa <- readRDS(pentads_file)
    } else {
        pentads_sa <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa") # HARD CODED
        saveRDS(pentads_sa, file.path(tempdir(), "sa_pentads.rds"))
    }

    # Add geometry
    summ_occu <- summ_occu %>%
        dplyr::left_join(pentads_sa %>%
                             dplyr::select(pentad), by = "pentad") %>%
        sf::st_sf()

    # Occupancy probabilities
    psi <- summ_occu %>%
        ggplot() +
        geom_sf(aes(fill = psi), color = NA) +
        scale_fill_viridis_c(limits = c(0, 1)) +
        ggtitle(paste(sp_name, year_sel)) +
        facet_wrap("lim")

    # Detection probabilities
    p <- summ_occu %>%
        dplyr::mutate(p = ifelse(is.na(p), 0, p)) %>%
        ggplot() +
        geom_sf(aes(fill = p), color = NA) +
        scale_fill_viridis_c(name = "p", limits = c(0, 1)) +
        ggtitle(paste(sp_name, year_sel)) +
        facet_wrap("lim")

    # Realized occupancy
    occu <- summ_occu %>%
        ggplot() +
        geom_sf(aes(fill = real_occu), color = NA) +
        scale_fill_viridis_c(limits = c(0, 1)) +
        ggtitle(paste(sp_name, year_sel)) +
        facet_wrap("lim")

    ggsave(file.path(config$out_dir, sp_code, paste0("occur_psi_", year_sel, "_", sp_code, ".png")), psi)
    ggsave(file.path(config$out_dir, sp_code, paste0("occur_p_", year_sel, "_", sp_code, ".png")), p)
    ggsave(file.path(config$out_dir, sp_code, paste0("occur_occu_", year_sel, "_", sp_code, ".png")), occu)

}