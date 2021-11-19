# In this script we fit an occupancy model with occuR for all species

library(occuR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(BIRDIE)

rm(list = ls())

birdie_dir <- "analysis/"

# For now, we want to select species present at Barberspan
bbpan <- BIRDIE::barberspan %>%
    pull(spp) %>%
    unique()

# Select a range of time. Occupancy models will be fitted from first to second
ini_year <- 2008
years <- c(ini_year, ini_year + 4)
years_ch <- paste(substring(as.character(years), 3, 4), collapse = "_")

for(i in seq_along(bbpan)){

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
    sp_detect <- SABAP::getSabapData(.spp_code = sp_sel,
                                     .region_type = "country",
                                     .region = "South Africa",
                                     .years = years[1]:years[2])


    # Load occupancy data -----------------------------------------------------

    # Load site data and subset years
    sitedata <- readRDS(paste0(birdie_dir, "data/site_dat_sa_gee_08_19.rds")) %>%
        dplyr::select(Name, lon, lat, watocc_ever, ends_with(match = as.character(years[1]:years[2])))

    # Load visit data, subset years and add detections
    visitdata <- readRDS(paste0(birdie_dir, "data/visit_dat_sa_gee_08_19.rds")) %>%
        filter(year %in% (years[1]:years[2])) %>%
        left_join(sp_detect %>%
                      dplyr::select(CardNo, obs = Spp) %>%
                      mutate(obs = if_else(obs == "-", 0, 1)),
                  by = "CardNo")


    # Format to occuR ---------------------------------------------------------

    occuRdata <- prepDataOccuR(sitedata, visitdata)

    # Remove data from missing pentads? (THIS SHOULDN'T HAPPEN WHEN PENTADS ARE DOWNLOADED FROM THE API)
    occuRdata$visit <- occuRdata$visit %>%
        dplyr::filter(!is.na(site))


    # Fit occupancy model -----------------------------------------------------

    # Define site model
    sitemod <- c("1", "s(watocc, bs = 'cs')", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')",
                 "t2(lon, lat, occasion, k = c(20, 3), bs = c('ts', 'cs'), d = c(2, 1))")

    # Define visit model
    visitmod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

    # Scale variables
    occuRdata$site <- occuRdata$site %>%
        mutate(tdiff = tmmx - tmmn,
               across(.col = c(prcp, tdiff), .fns = ~scale(.x)))

    print(paste0("Fitting model at ", Sys.time(), ". This will take a while..."))

    # Smooth for spatial effect on psi
    fit <- fit_occu(forms = list(reformulate(visitmod, response = "p"),
                                 reformulate(sitemod, response = "psi")),
                    visit_data = occuRdata$visit,
                    site_data = occuRdata$site,
                    print = TRUE)

    saveRDS(fit, paste0("analysis/output/occur_fit_", years_ch, "_", sp_sel, ".rds"))


    # Predict occupancy -------------------------------------------------------

    # Extract scaling factors
    sc <- lapply(occuRdata$site, attributes)
    sc <- sc[!sapply(sc, is.null)]

    # Prepare data to predict psi and p
    gm <- sitedata %>%
        dplyr::select(Name)

    pred_data <- sitedata %>%
        sf::st_drop_geometry()  %>%
        group_by(Name) %>%
        mutate(site = cur_group_id()) %>%
        ungroup() %>%
        pivot_longer(cols = -c(id, Name, site, lon, lat, water)) %>%
        tidyr::separate(name, into = c("covt", "year"), sep = "_") %>%
        pivot_wider(names_from = covt, values_from = value) %>%
        mutate(year = as.numeric(year),
               tdiff = tmax - tmin) %>%
        dplyr::filter(year >= years[1], year <= years[2]) %>%
        mutate(prcp = scale(prcp,
                            center = sc$prcp$`scaled:center`,
                            scale = sc$prcp$`scaled:scale`),
               tdiff = scale(tdiff,
                             center = sc$tdiff$`scaled:center`,
                             scale = sc$tdiff$`scaled:scale`)) %>%
        group_by(year) %>%
        mutate(occasion = cur_group_id()) %>%
        ungroup() %>%
        data.table::as.data.table()

    # Predict
    pred <- predict(fit, occuRdata$visit,  pred_data, nboot = 1000)


    # Estimate realized occupancy ---------------------------------------------

    # Calculate probability of non-detections for each pentad visited
    p_nondet <- occuRdata$visit %>%
        dplyr::select(Pentad, site, occasion, obs) %>%
        mutate(ub = apply(pred$pboot, 2, quantile, 0.975),
               lb = apply(pred$pboot, 2, quantile, 0.025),
               med = apply(pred$pboot, 2, quantile, 0.5),
               est = pred$p) %>%
        pivot_longer(cols = -c(Pentad, site, occasion, obs),
                     names_to = "lim", values_to = "p") %>%
        group_by(Pentad, site, occasion, lim) %>%
        summarize(pp = prod(1-p),
                  obs = max(obs))

    # From probability of non-detection calculate the conditional occupancy probs
    # and plot
    pred_occu <- pred_data %>%
        as.data.frame() %>%
        left_join(gm, by = "Name") %>%
        sf::st_sf() %>%
        mutate(ub = apply(pred$psiboot, 2, quantile, 0.975),
               lb = apply(pred$psiboot, 2, quantile, 0.025),
               med = apply(pred$psiboot, 2, quantile, 0.5),
               est = pred$psi[,1]) %>%
        pivot_longer(cols = c("ub", "lb", "med", "est"),
                     names_to = "lim", values_to = "psi") %>%
        left_join(p_nondet, by = c("Name" = "Pentad", "occasion", "lim")) %>%
        mutate(real_occu = case_when(obs == 1 ~ 1,
                                     is.na(obs) ~ psi,
                                     obs == 0 ~ psi*pp / (1 - psi + psi*pp))) %>%
        dplyr::select(Name, year, psi, pp, real_occu, lim)

    for(t in seq(years[1], years[2], 1)){
        # save data if year < 2010 or otherwise if the year is in the middle
        # of the series or higher (middle should give the most accurate temporal
        # estimate)
        if(t < 2010 | (years[2] - t) < 2){
            year_sel <- substring(as.character(t), 3, 4)
            pred_occu %>%
                sf::st_drop_geometry() %>%
                dplyr::filter(year == t) %>%
                write.csv(paste0("analysis/output/occur_pred_", year_sel, "_", sp_sel, ".csv"),
                          row.names = FALSE)
        }
    }


    # Plot model -----------------------------------------------------------

    for(t in seq(years[1], years[2], 1)){
        # save data if year < 2010 or otherwise if the year is in the middle
        # of the series or higher (middle should give the most accurate temporal
        # estimate)
        if(t < 2010 | (years[2] - t) < 2){
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
                mutate(pp = if_else(is.na(pp), 1, pp)) %>%
                ggplot() +
                geom_sf(aes(fill = 1 - pp), size = 0.01) +
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

            ggsave(paste0("analysis/output/occur_psi_", year_sel, "_", sp_sel, ".png"), psi)
            ggsave(paste0("analysis/output/occur_p_", year_sel, "_", sp_sel, ".png"), p)
            ggsave(paste0("analysis/output/occur_occu_", year_sel, "_", sp_sel, ".png"), occu)
        }
    }
}
