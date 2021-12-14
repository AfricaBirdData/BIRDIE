# In this script we fit an occupancy model with occuR for all species

library(occuR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(BIRDIE)

rm(list = ls())

# Define data and output directories
data_dir <- "analysis/data/"
out_dir <- "/drv_birdie/birdie_ftp/"

# For now, we want to select species present at Barberspan
bbpan <- BIRDIE::barberspan %>%
    pull(spp) %>%
    unique()

# Select a range of time. Occupancy models will be fitted from first to second
ini_year <- 2008
year_range <- c(ini_year, ini_year + 4)
years_ch <- paste(substring(as.character(year_range), 3, 4), collapse = "_")
years <- year_range[1]:year_range[2]

# Load site data
sitedata <- readRDS(paste0(data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    dplyr::select(Name, lon, lat, watocc_ever, dist_coast, ends_with(match = as.character(year_sel)))

# I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE
sitedata <- sitedata %>%
    tidyr::drop_na()


# Define models -----------------------------------------------------------

# Detection
visit_mod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

# Occupancy

# For inland species
land_mod <- c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')",
              "t2(lon, lat, occasion, k = c(20, 4), bs = c('ts', 'cs'), d = c(2, 1))")

# For coastal species
coast_mod <- c("1", "dist_coast", "prcp", "tdiff", "watext", "watrec",
               "t2(lon, lat, occasion, k = c(20, 4), bs = c('ts', 'cs'), d = c(2, 1))")


for(i in seq_along(bbpan)){

    # i = 92 this is a coastal species

    # Select a species and a region -------------------------------------------

    sp_sel <- bbpan[i]

    print(paste0("Working on species ", sp_sel, " (", i, " of ", length(bbpan), ")"))

    sp_name <- BIRDIE::barberspan %>%
        dplyr::filter(spp == sp_sel) %>%
        mutate(name = paste(taxon.Common_species, taxon.Common_group)) %>%
        pull(name) %>%
        unique()

    # Download species detection
    print("Downloading from SABAP")
    sp_detect <- ABAP::getAbapData(.spp_code = sp_sel,
                                   .region_type = "country",
                                   .region = "South Africa",
                                   .years = years)


    # Load occupancy data -----------------------------------------------------

    # Load visit data, subset years and add detections
    visitdata <- readRDS(paste0(data_dir, "visit_dat_sa_gee_08_19.rds")) %>%
        filter(year %in% years) %>%
        left_join(sp_detect %>%
                      dplyr::select(CardNo, obs = Spp) %>%
                      mutate(obs = if_else(obs == "-", 0, 1)),
                  by = "CardNo")


    # Format to occuR ---------------------------------------------------------

    occuRdata <- prepDataOccuR(site_data = sitedata %>%
                                   sf::st_drop_geometry() %>%
                                   gatherYearFromVars(vars = names(.)[-c(1:6)], sep = "_") %>%
                                   mutate(tdiff = tmmx - tmmn),
                               visit_data = visitdata,
                               scaling = list(visit = NULL,
                                              site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")))

    # Remove data from missing pentads? (THIS SHOULDN'T HAPPEN WHEN PENTADS ARE DOWNLOADED FROM THE API)
    occuRdata$visit <- occuRdata$visit %>%
        dplyr::filter(!is.na(site))

    occuRdata$site <- occuRdata$site %>%
        dplyr::filter(site %in% unique(occuRdata$visit$site))


    # Fit occupancy model -----------------------------------------------------

    # Determine if the species is coastal
    coast <- isSpCoastal(sp_sel, out_dir, reformulate(visit_mod, response = "p"))

    if(coast){
        site_mod <- coast_mod
    } else {
        site_mod <- land_mod
    }

    print(paste0("Fitting model at ", Sys.time(), ". This will take a while..."))

    # Smooth for spatial effect on psi
    tryCatch({
        fit <- fit_occu(forms = list(reformulate(visit_mod, response = "p"),
                                     reformulate(site_mod, response = "psi")),
                        visit_data = occuRdata$visit,
                        site_data = occuRdata$site,
                        print = TRUE)

        saveRDS(fit, paste0(out_dir, sp_sel, "/occur_fit_", years_ch, "_", sp_sel, ".rds"))
    }, error = function(e){
        sink(paste0(out_dir, sp_sel, "/failed_fit_", sp_sel,".txt"))
        print(e)
        sink()}) # TryCatch fit


    # Predict occupancy -------------------------------------------------------
    tryCatch({
        pred_occu <- predictOccuR(fit, occuRdata, sitedata,
                                  years = years,
                                  scaling = TRUE)
        # Save predictions
        for(t in years){
            # save data if year < 2010 or otherwise if the year is in the middle
            # of the series or higher (middle should give the most accurate temporal
            # estimate)
            if(t < 2010 | (year_range[2] - t) < 2){
                year_sel <- substring(as.character(t), 3, 4)
                pred_occu %>%
                    sf::st_drop_geometry() %>%
                    dplyr::filter(year == t) %>%
                    write.csv(paste0(out_dir, sp_sel, "/occur_pred_", year_sel, "_", sp_sel, ".csv"),
                              row.names = FALSE)
            }
        }


        # Plot model -----------------------------------------------------------

        for(t in years){
            # save data if year < 2010 or otherwise if the year is in the middle
            # of the series or higher (middle should give the most accurate temporal
            # estimate)
            if(t < 2010 | (year_range[2] - t) < 2){
                year_sel <- substring(as.character(t), 3, 4)

                # Occupancy probabilities
                psi <- pred_occu %>%
                    dplyr::filter(year == t) %>%
                    ggplot() +
                    geom_sf(aes(fill = psi), size = 0.01) +
                    scale_fill_viridis_c(limits = c(0, 1)) +
                    ggtitle(sp_name) +
                    facet_wrap("lim")

                # Detection probabilities
                p <- pred_occu %>%
                    dplyr::filter(year == t) %>%
                    mutate(p = if_else(is.na(p), 0, p)) %>%
                    ggplot() +
                    geom_sf(aes(fill = p), size = 0.01) +
                    scale_fill_viridis_c(name = "p", limits = c(0, 1)) +
                    ggtitle(sp_name) +
                    facet_wrap("lim")

                # Realized occupancy
                occu <- pred_occu %>%
                    dplyr::filter(year == t) %>%
                    ggplot() +
                    geom_sf(aes(fill = real_occu), size = 0.01) +
                    scale_fill_viridis_c(limits = c(0, 1)) +
                    ggtitle(sp_name) +
                    facet_wrap("lim")

                ggsave(paste0(out_dir, sp_sel, "/occur_psi_", year_sel, "_", sp_sel, ".png"), psi)
                ggsave(paste0(out_dir, sp_sel, "/occur_p_", year_sel, "_", sp_sel, ".png"), p)
                ggsave(paste0(out_dir, sp_sel, "/occur_occu_", year_sel, "_", sp_sel, ".png"), occu)
            }
        }
    }, error = function(e){
        sink(paste0(out_dir, sp_sel, "/failed_pred_", sp_sel,".txt"))
        print(e)
        sink()}) # TryCatch predict


}
