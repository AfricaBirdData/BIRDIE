#' Estimate occupancy for species type
#'
#' @inheritParams ppl_estimate_distr_sp_type
#'
#' @return
#' @import sf
#'
#' @examples
estimateOccuSpType <- function(sp_type, year, config){

    # Lookup species in the group
    if(sp_type == "group"){
        species <- config$species[1:4] # This needs to be fixed when groups are properly defined
    } else {
        stop("Sorry the group doesn't exist")
    }

    # Save predictions
    for(t in seq_along(config$years)){

        # save data and plots if the year is in the middle of the series or
        # higher (middle should give the most accurate temporal estimate)

        if((t > config$dyear/2) | (config$year < (2009 + config$dyear/2))){

            yy <- substring(as.character(config$years[t]), 3, 4)

            # Retrieve occupancy prediction files
            ffs <- list()

            for(i in seq_along(species)){
                ffs[[i]] <- read.csv(file.path(config$fit_dir, species[i], paste0("occur_pred_", yy, "_", species[i], ".csv")))
            }

            ffs <- bind_rows(ffs)

            ffs <- ffs %>%
                dplyr::mutate(q = 1 - real_occu) %>%
                dplyr::group_by(Pentad, lim) %>%
                dplyr::summarise(q = prod(q)) %>%
                dplyr::mutate(real_occu = 1 - q) %>%
                dplyr::ungroup()

            ## SAVE
            ffs %>%
                write.csv(file.path(config$fit_dir, sp_type, paste0("occur_pred_", yy, "_", sp_type, ".csv")),
                          row.names = FALSE)

            ## PLOT

            # Add geometry
            gm <- readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")) %>%
                dplyr::select(Pentad = Name)

            pred_sel <- gm %>%
                dplyr::left_join(ffs, by = "Pentad") %>%
                dplyr::filter(!is.na(real_occu))


            # Realized occupancy
            occu <- pred_sel %>%
                ggplot() +
                geom_sf(aes(fill = real_occu), size = 0.01) +
                scale_fill_viridis_c(limits = c(0, 1)) +
                ggtitle(paste(sp_type, config$years[t])) +
                facet_wrap("lim")

            ggsave(file.path(config$fit_dir, sp_type, paste0("occur_occu_", yy, "_", sp_type, ".png")), occu)

        }
    }
}
