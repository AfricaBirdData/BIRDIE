# In this script we fit an occupancy model with occuR for all species

library(occuR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(BIRDIE)
library(sf)

rm(list = ls())


# Script parameters -------------------------------------------------------

server <- FALSE
ini_year <- 2008

if(server){

    # Define data and output directories
    data_dir <- "/home/birdie/analysis/data/"
    out_dir <- "/drv_birdie/birdie_ftp/"

    # Define years to fit
    dyear <- 4  # this will be 4 in the server and 3 locally
    year_range <- c(ini_year, ini_year + dyear)

    # and spatio-temporal effect
    sptemp <- "t2(lon, lat, occasion, k = c(20, 4), bs = c('ts', 'cs'), d = c(2, 1))" # this will be c(20, 4) in the server and c(15, 3) locally

    # Define species to fit models to
    species <- BIRDIE::barberspan %>% # For now, we want to select species present at Barberspan
        pull(SppRef) %>%
        unique()


} else {

    # Define data and output directories
    data_dir <- "analysis/data/"
    out_dir <- "analysis/out_nosync/"

    # Define years to fit
    dyear <- 2  # this will be 4 in the server and 2 locally
    year_range <- c(ini_year, ini_year + dyear)

    # and spatio-temporal effect
    sptemp <- "t2(lon, lat, occasion, k = c(15, 3), bs = c('ts', 'cs'), d = c(2, 1))" # this will be c(20, 4) in the server and c(15, 3) locally

    # Define species to fit models to
    species <- c(4, 6, 41, 235, 240)

}

# Define a range of years to fit occupancy model
years_ch <- paste(substring(as.character(year_range), 3, 4), collapse = "_")
years <- year_range[1]:year_range[2]


# Load site data ----------------------------------------------------------

sitedata <- readRDS(paste0(data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    dplyr::select(Name, lon, lat, watocc_ever, dist_coast, ends_with(match = as.character(years)))

# I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE
sitedata <- sitedata %>%
    tidyr::drop_na()


# Define models -----------------------------------------------------------

# Detection
visit_mod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

# Occupancy
site_mods <- list(mod1 = c("-1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec", sptemp),
                  mod2 = c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"),
                  mod3 = c("1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec"))



# Predict from models -----------------------------------------------------

for(i in seq_along(species)){

    # This will be used for memory cleaning
    keep <- ls()

    # Select a species and a region -------------------------------------------

    sp_sel <- species[i]

    print(paste0("Working on species ", sp_sel, " (", i, " of ", length(species), ")"))

    sp_name <- BIRDIE::barberspan %>%
        dplyr::filter(SppRef == sp_sel) %>%
        mutate(name = paste(Common_species, Common_group)) %>%
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
                                   gatherYearFromVars(vars = names(.)[-c(1:5)], sep = "_") %>% # check that 1:5 are the variables that don't change over time
                                   mutate(tdiff = tmmx - tmmn),
                               visit_data = visitdata,
                               scaling = list(visit = NULL,
                                              site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")))

    # Remove data from missing pentads? (THIS SHOULDN'T HAPPEN WHEN PENTADS ARE DOWNLOADED FROM THE API)
    occuRdata$visit <- occuRdata$visit %>%
        dplyr::filter(!is.na(site))

    occuRdata$site <- occuRdata$site %>%
        dplyr::filter(site %in% unique(occuRdata$visit$site))


    # Predict occupancy -------------------------------------------------------

    # Load fit
    fit <- readRDS(paste0(out_dir, sp_sel, "/occur_fit_", years_ch, "_", sp_sel, ".rds"))

    # Prepare prediction data
    pred_data <- prepPredictDataOccuR(occuRdata, sf::st_drop_geometry(sitedata),
                                      years = years, scaling = TRUE)

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
        for(t in seq_along(years)){

            # save data and plots if the year is in the middle of the series or
            # higher (middle should give the most accurate temporalestimate)

            if(t > dyear/2){

                # select year
                year_sel <- substring(as.character(years[t]), 3, 4)

                # Subset predictions
                pred_sel <- summ_occu %>%
                    dplyr::filter(year == years[t])

                # Save predictions
                pred_sel %>%
                    write.csv(paste0(out_dir, sp_sel, "/occur_pred_", year_sel, "_", sp_sel, ".csv"),
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
                    ggtitle(paste(sp_name, years[t])) +
                    facet_wrap("lim")

                # Detection probabilities
                p <- pred_sel %>%
                    mutate(p = if_else(is.na(p), 0, p)) %>%
                    ggplot() +
                    geom_sf(aes(fill = p), size = 0.01) +
                    scale_fill_viridis_c(name = "p", limits = c(0, 1)) +
                    ggtitle(paste(sp_name, years[t])) +
                    facet_wrap("lim")

                # Realized occupancy
                occu <- pred_sel %>%
                    ggplot() +
                    geom_sf(aes(fill = real_occu), size = 0.01) +
                    scale_fill_viridis_c(limits = c(0, 1)) +
                    ggtitle(paste(sp_name, years[t])) +
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

    # Clean
    rm(list = setdiff(ls(), keep))
    gc()

}
