library(occuR)
library(tidyverse)
library(BIRDIE)

rm(list = ls())


# For now, we want to select species present at Barberspan
bbpan <- BIRDIE::barberspan %>%
    pull(spp) %>%
    unique()


# Select a species and a region -------------------------------------------

# We are most interested in the Maccoa Duck VU (103) and the Cape Cormorant EN (48)
i <- 1
sp_sel <- bbpan[i]
sp_sel <- 103

sp_name <- BIRDIE::barberspan %>%
    dplyr::filter(spp == sp_sel) %>%
    mutate(name = paste(taxon.Common_species, taxon.Common_group)) %>%
    pull(name) %>%
    unique()


# Load occupancy data -----------------------------------------------------

# Load site data
sitedata <- readRDS("analysis/data/site_dat_sa_wcovts_16_19.rds")

# Load visit data
visitdata <- readRDS(paste0("analysis/data/visit_dat_", sp_sel, "_wcovts_16_19.rds"))


# Format to occuR ---------------------------------------------------------

occuRdata <- prepDataOccuR(sitedata, visitdata)


# Predict occupancy -------------------------------------------------------

fit <- readRDS(paste0("analysis/output/", sp_sel, "_occur_fit_16_19.rds"))

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
    mutate(year = as.numeric(year)) %>%
    group_by(year) %>%
    mutate(occasion = cur_group_id()) %>%
    ungroup() %>%
    data.table::as.data.table()

# Predict
pred <- predict(fit, occuRdata$visit,  pred_data, nboot = 1000)


# Estimate realized occupancy ---------------------------------------------

# Calculate probability of non-detections for each pentad visited
p_nondet <- occuRdata$visit %>%
    mutate(p = pred$p) %>%
    group_by(Pentad, site, occasion) %>%
    summarize(pp = prod(1-p),
              obs = max(obs))

# From probability of non-detection calculate the conditional occupancy probs
# and plot
pred_occu <- pred_data %>%
    as.data.frame() %>%
    left_join(gm, by = "Name") %>%
    sf::st_sf() %>%
    mutate(psi = pred$psi) %>%
    left_join(p_nondet, by = c("Name" = "Pentad", "occasion")) %>%
    mutate(real_occu = case_when(obs == 1 ~ 1,
                                 is.na(obs) ~ psi,
                                 obs == 0 ~ psi*pp / (1 - psi + psi*pp))) %>%
    dplyr::select(Name, year, psi, pp, real_occu)


# Explore model -----------------------------------------------------------

# Occupancy probabilities
pred_occu %>%
    ggplot() +
    geom_sf(aes(fill = psi), size = 0.01) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    ggtitle(sp_name) +
    facet_wrap("year")

# Detection probabilities
pred_occu %>%
    mutate(pp = if_else(is.na(pp), 1, pp)) %>%
    ggplot() +
    geom_sf(aes(fill = 1 - pp), size = 0.01) +
    scale_fill_viridis_c(name = "p", limits = c(0, 1)) +
    ggtitle(sp_name) +
    facet_wrap("year")

# Realized occupancy
pred_occu %>%
    ggplot() +
    geom_sf(aes(fill = real_occu), size = 0.01) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    ggtitle(sp_name) +
    facet_wrap("year")

# ggsave(filename = paste0("analysis/out_nosync/occu_", sp_sel, "_15_19.png"))

sitedata %>%
    ggplot() +
    geom_sf(aes(fill = log(water+0.1)), size = 0.01) +
    scale_fill_viridis_c()

pred_data %>%
    as.data.frame() %>%
    left_join(gm, by = "Name") %>%
    sf::st_sf() %>%
    ggplot() +
    geom_sf(aes(fill = prcp), size = 0.01) +
    scale_fill_viridis_c() +
    facet_wrap("year")

# Plot detection
det <- occuRdata$visit %>%
    group_by(Pentad, year) %>%
    summarise(dets = sum(obs)) %>%
    ungroup() %>%
    mutate(det = if_else(dets > 0, 1L, 0L))

gm %>%
    right_join(det, by = c("Name" = "Pentad")) %>%
    ggplot() +
    geom_sf(aes(fill = factor(det)), size = 0.01) +
    scale_fill_viridis_d(option = "E") +
    facet_wrap("year")


# Explore non-linear effects ----------------------------------------------

pred_data <- pred_data %>%
    mutate(prcp = 0,
           watext = 0,
           watrec = 0,
           dist_coast = 0,
           ndvi = 0,
           tmax = 0,
           tmin = 0) %>%
    data.table::as.data.table()

pred <- predict(fit, occuRdata$visit,  pred_data, nboot = 0)

pred_data %>%
    as.data.frame() %>%
    left_join(gm, by = "Name") %>%
    sf::st_sf() %>%
    mutate(psi = pred$psi) %>%
    ggplot() +
    geom_sf(aes(fill = psi), size = 0.01) +
    scale_fill_viridis_c() +
    ggtitle(sp_name) +
    facet_wrap("year")

pred_data %>%
    as.data.frame() %>%
    left_join(gm, by = "Name") %>%
    sf::st_sf() %>%
    mutate(psi = log(pred$psi/(1-pred$psi)) - fit$fit$par["beta_psi"]) %>%
    ggplot() +
    geom_sf(aes(fill = psi), size = 0.01) +
    scale_fill_viridis_c() +
    facet_wrap("year")



