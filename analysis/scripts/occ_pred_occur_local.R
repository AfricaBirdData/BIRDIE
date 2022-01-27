# In this script we fit an occupancy model with occuR for all species

library(occuR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(BIRDIE)
library(sf)

rm(list = ls())


# Script parameters -------------------------------------------------------

year_sel <- 2009

config <- configPreambOccuR(year = year_sel, server = FALSE)


# Load site and visit data ------------------------------------------------

# Load data and subset years
sitedata <- readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    dplyr::select(Name, lon, lat, watocc_ever, dist_coast, ends_with(match = as.character(config$years))) %>%
    tidyr::drop_na()   # I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE

visitdata <- readRDS(file.path(config$data_dir, "visit_dat_sa_gee_08_19.rds")) %>%
    filter(year %in% config$years)


# Predict from models -----------------------------------------------------

for(i in seq_along(config$species)){

    # This will be used for memory cleaning
    keep <- ls()

    # Select a species and a region -------------------------------------------

    sp_sel <- config$species[i]

    print(paste0("Working on species ", sp_sel, " (", i, " of ", length(config$species), ")"))

    sp_name <- BIRDIE::barberspan %>%
        dplyr::filter(SppRef == sp_sel) %>%
        mutate(name = paste(Common_species, Common_group)) %>%
        pull(name) %>%
        unique()


    # Format to occuR ---------------------------------------------------------

    occuRdata <- BIRDIE::prepDataOccuR(spp_code = sp_sel,
                                       years = config$years,
                                       site_data = sitedata,
                                       visit_data = visitdata,
                                       download = TRUE)


    # Predict occupancy -------------------------------------------------------

    # Load fit
    fit <- with(config,
                readRDS(file.path(config$fit_dir, sp_sel, paste0("occur_fit_", years_ch, "_", sp_sel, ".rds"))))

    # Prepare prediction data
    pred_data <- prepPredictDataOccuR(occuRdata, sf::st_drop_geometry(sitedata),
                                      years = config$years, scaling = TRUE)

    print(paste0("Predicting from model at ", Sys.time()))

    tryCatch({

        # Predict
        pred_occu <- predict(fit, occuRdata$visit, pred_data, nboot = 1000)

        # Summarize predictions
        summ_occu <- summarizePredOccuR(pred_p = pred_occu$pboot,
                                        pred_psi = pred_occu$psiboot,
                                        pred_data = pred_data,
                                        visit_data = occuRdata$visit,
                                        quants = c(0.025, 0.5, 0.975))

        # Save predictions
        for(t in seq_along(config$years)){

            # save data and plots if the year is in the middle of the series or
            # higher (middle should give the most accurate temporal estimate)

            if((t > config$dyear/2) | (config$year < (2009 + config$dyear/2))){

                # select year
                yy <- substring(as.character(config$years[t]), 3, 4)

                # Subset predictions and add species
                pred_sel <- summ_occu %>%
                    dplyr::filter(year == config$years[t]) %>%
                    dplyr::mutate(species = sp_sel) %>%
                    dplyr::select(species, everything())

                # Save predictions
                pred_sel %>%
                    write.csv(file.path(config$fit_dir, sp_sel, paste0("occur_pred_", yy, "_", sp_sel, ".csv")),
                              row.names = FALSE)

                ## PLOTS

                # Add geometry
                pred_sel <- sitedata %>%
                    dplyr::select(Name) %>%
                    dplyr::left_join(pred_sel, by = "Name")

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

                ggsave(file.path(config$fit_dir, sp_sel, paste0("occur_psi_", yy, "_", sp_sel, ".png")), psi)
                ggsave(file.path(config$fit_dir, sp_sel, paste0("occur_p_", yy, "_", sp_sel, ".png")), p)
                ggsave(file.path(config$fit_dir, sp_sel, paste0("occur_occu_", yy, "_", sp_sel, ".png")), occu)

            }
        }

    }, error = function(e){
        sink(file.path(config$fit_dir, sp_sel, paste0("failed_pred_", sp_sel,".txt")))
        print(e)
        sink()}) # TryCatch predict

    # Clean
    rm(list = setdiff(ls(), keep))
    gc()

}
