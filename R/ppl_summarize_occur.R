#' Summarize prediction from occuR model
#'
#' @inheritParams ppl_run_pipe_distr
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
ppl_summarize_occur <- function(sp_code, sp_name, year, config, ...){

    # Predict from model
    pred_occu <- predictOccuR(sp_code, year, config, ...)

    print("Summarizing predictions")

    # Summarize predictions
    summ_occu <- summarizeOccuRPred(pred_p = pred_occu$pred$pboot,
                                    pred_psi = pred_occu$pred$psiboot,
                                    pred_data = pred_occu$pred_data,
                                    visit_data = pred_occu$occuRdata$visit,
                                    quants = c(0.025, 0.5, 0.975))

    # Save predictions
    for(t in seq_along(config$years)){

        # save data and plots if the year is in the middle of the series or
        # higher (middle should give the most accurate temporal estimate)

        if((t > config$dur/2) | (config$year < (2009 + config$dur))){

            # select year
            yy <- substring(as.character(config$years[t]), 3, 4)

            # Subset predictions and add species
            pred_sel <- summ_occu %>%
                dplyr::filter(year == config$years[t]) %>%
                dplyr::mutate(species = sp_code) %>%
                dplyr::select(species, everything())

            # Save predictions
            pred_sel %>%
                write.csv(file.path(config$out_dir, sp_code, paste0("occur_pred_", yy, "_", sp_code, ".csv")),
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
            pred_sel <- pred_occu$pred_sites %>%
                dplyr::left_join(pentads_sa, by = c("Pentad" = "pentad")) %>%
                sf::st_as_sf() %>%
                dplyr::select(Pentad) %>%
                dplyr::left_join(pred_sel, by = "Pentad")

            # Occupancy probabilities
            psi <- pred_sel %>%
                ggplot() +
                geom_sf(aes(fill = psi), size = 0.01) +
                scale_fill_viridis_c(limits = c(0, 1)) +
                ggtitle(paste(sp_name, config$years[t])) +
                facet_wrap("lim")

            # Detection probabilities
            p <- pred_sel %>%
                mutate(p = if_else(is.na(p), 0, p)) %>%
                ggplot() +
                geom_sf(aes(fill = p), size = 0.01) +
                scale_fill_viridis_c(name = "p", limits = c(0, 1)) +
                ggtitle(paste(sp_name, config$years[t])) +
                facet_wrap("lim")

            # Realized occupancy
            occu <- pred_sel %>%
                ggplot() +
                geom_sf(aes(fill = real_occu), size = 0.01) +
                scale_fill_viridis_c(limits = c(0, 1)) +
                ggtitle(paste(sp_name, config$years[t])) +
                facet_wrap("lim")

            ggsave(file.path(config$out_dir, sp_code, paste0("occur_psi_", yy, "_", sp_code, ".png")), psi)
            ggsave(file.path(config$out_dir, sp_code, paste0("occur_p_", yy, "_", sp_code, ".png")), p)
            ggsave(file.path(config$out_dir, sp_code, paste0("occur_occu_", yy, "_", sp_code, ".png")), occu)

        }
    }

}
