library(occuR)
library(tidyverse)
library(BIRDIE)

rm(list = ls())


# For now, we want to select species present at Barberspan
bbpan <- BIRDIE::barberspan %>%
    pull(spp) %>%
    unique()


# Select a species and a region -------------------------------------------

i <- 1

sp_sel <- bbpan[i]

sp_name <- BIRDIE::barberspan %>%
    dplyr::filter(spp == sp_sel) %>%
    mutate(name = paste(taxon.Common_species, taxon.Common_group)) %>%
    pull(name) %>%
    unique()


# Load occupancy site data -------------------------------------------------

# Load site data
sitedata <- readRDS("analysis/data/site_dat_sa_wcovts.rds")


# Prepare occupancy visit data --------------------------------------------

future::plan("multisession", workers = 6)

visitdata <- prepOccVisitData(region = "South Africa",
                              sites = sitedata,
                              species = sp_sel,
                              years = 2008:2011,
                              clim_covts = c("prcp", "tmax", "tmin", "aet", "pet"),
                              covts_dir = "analysis/downloads/",
                              file_fix = c("terraClim_", "_03_19"),
                              savedir = "analysis/data/pentads_sa.rds")

future::plan("sequential")


# Format to occuR ---------------------------------------------------------

occuRdata <- prepDataOccuR(sitedata, visitdata)


# Fit occupancy model -----------------------------------------------------

# Define site model
sitemod <- c("1", "s(water, bs = 'cs')", "s(prcp, bs = 'cs')", "s(tmax - tmin, bs = 'cs')",
             "t2(lon, lat, occasion, k = c(15, 3), bs = c('ts', 'cs'), d = c(2,1))")
# The original objective was 25 knots for the spatial effect (although I ran into memory issues)

# Define visit model
visitmod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

# Scale variables
# visitvars <- visitvars %>%
#     mutate(across(.col = -c(lon, lat, year, month, Pentad, obs, site, occasion, visit), .fns = ~scale(.x)))

# Smooth for spatial effect on psi
fit <- fit_occu(forms = list(reformulate(visitmod, response = "p"),
                             reformulate(sitemod, response = "psi")),
                visit_data = occuRdata$visit,
                site_data = occuRdata$site)

saveRDS(fit, paste0("analysis/output/", sp_sel, "_occur_fit.rds"))



# Predict occupancy -------------------------------------------------------

fit <- readRDS(paste0("analysis/output/", sp_sel, "_occur_fit.rds"))

# Predict psi
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
    dplyr::select(Name, year, psi, real_occu)


# Explore model -----------------------------------------------------------

# Occupancy probabilities
pred_occu %>%
    ggplot() +
    geom_sf(aes(fill = psi), size = 0.01) +
    scale_fill_viridis_c() +
    ggtitle(sp_name) +
    facet_wrap("year")

# Realized occupancy
pred_occu %>%
    ggplot() +
    geom_sf(aes(fill = real_occu), size = 0.01) +
    scale_fill_viridis_c() +
    ggtitle(sp_name) +
    facet_wrap("year")

# ggsave(filename = "analysis/out_nosync/occu_6_08_11.png")

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
    mutate(prcp = min(prcp),
           water = min(water),
           tmax = min(tmax),
           tmin = min(tmin)) %>%
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
    mutate(psi = log(pred$psi/(1-pred$psi)) - 2.0150) %>%
    ggplot() +
    geom_sf(aes(fill = psi), size = 0.01) +
    scale_fill_viridis_c() +
    facet_wrap("year")



