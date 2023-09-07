#' Create prediction object directly from ABAP data
#'
#' @description
#' This function is used when it is not possible to run an occupancy model for
#' an species and year. Then, the raw data is presented but they must have the
#' same structure as the prediction object to be able to integrate them in the
#' database.
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
createPredFromAbap <- function(sp_code, year_sel, config){

    # Download geometry if not present on disk
    pentads_file <- file.path(tempdir(), "sa_pentads.rds")

    if(file.exists(pentads_file)){
        pentads_sa <- readRDS(pentads_file)
    } else {
        pentads_sa <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa") # HARD CODED
        saveRDS(pentads_sa, file.path(tempdir(), "sa_pentads.rds"))
    }

    # Load detection data
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))
    det_data <- utils::read.csv(detfile, check.names = FALSE)

    # Summarise detection data
    det_data <- det_data %>%
        dplyr::group_by(Pentad) %>%
        dplyr::summarise(dets = sum(obs)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(dets = ifelse(dets == 0, 0, 1))

    # Join pentads and detections
    det_pentad <- pentads_sa %>%
        dplyr::left_join(det_data, by = c("pentad" = "Pentad"))

    # Form dataframe with same characteristics as predictions
    out <- det_pentad %>%
        # sf::st_drop_geometry() %>%
        dplyr::mutate(species = sp_code,
                      year = year_sel,
                      psi = NA,
                      p = NA,
                      real_occu_lb = NA,
                      real_occu_med = NA,
                      real_occu_ub = NA,
                      real_occu_est = dets) %>%
        tidyr::pivot_longer(col = c("real_occu_lb", "real_occu_med", "real_occu_ub", "real_occu_est"),
                            names_to = "lim", values_to = "real_occu") %>%
        dplyr::select(species, pentad, year, psi, p, real_occu, lim) %>%
        dplyr::mutate(lim = gsub("real_occu_", "", lim))

    # Plot

    # Species name for plot titles
    sp_name <- BIRDIE::waterbirds %>%
        dplyr::filter(SppRef == sp_code) %>%
        dplyr::mutate(name = paste(Common_species, Common_group)) %>%
        dplyr::mutate(name = gsub(" NA|NA ", "", name)) %>% # in case there are NAs in species or group
        dplyr::pull(name) %>%
        unique()

    occu_plot <- out %>%
        dplyr::filter(lim == "est") %>%
        ggplot() +
        geom_sf(aes(fill = real_occu), color = NA) +
        scale_fill_viridis_c(limits = c(0, 1)) +
        ggtitle(paste(sp_name, year_sel)) +
        facet_wrap("lim")

    return(list(df = out %>%
                    sf::st_drop_geometry(),
                occu_plot = occu_plot))


}
