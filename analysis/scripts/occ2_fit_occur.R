library(occuR)
library(tidyverse)

rm(list = ls())


# Select a species and a region -------------------------------------------

sp_sel <- 6


# Load occupancy data -----------------------------------------------------

# Load site data
sitedata <- readRDS("analysis/data/site_dat_sa_wcovts.rds")

# Load visit data
visitdata <- readRDS(paste0("analysis/data/visit_dat_", sp_sel, "_wcovts.rds"))


# Format to occuR ---------------------------------------------------------

# Visited sites
visits <- visitdata %>%
    distinct(SiteName, year) %>%
    mutate(keep = 1)

sites <- unique(visitdata$SiteName)

site_data <- sitedata %>%
    dplyr::select(-id) %>%
    sf::st_drop_geometry()  %>%
    tidyr::pivot_longer(cols = -c(Name, lon, lat, water)) %>%
    tidyr::separate(name, into = c("covt", "year"), sep = "_") %>%
    pivot_wider(names_from = covt, values_from = value) %>%
    mutate(year = as.numeric(year)) %>%
    left_join(visits, by = c("Name" = "SiteName", "year" = "year")) %>%
    filter(keep == 1) %>%
    dplyr::select(-keep) %>%
    arrange(lat, lon) %>%
    group_by(Name) %>%
    mutate(site = cur_group_id()) %>%
    ungroup() %>%
    group_by(year) %>%
    mutate(occasion = cur_group_id()) %>%
    ungroup() %>%
    data.table::as.data.table()

visit_data <- visitdata %>%
    filter(year %in% unique(site_data$year)) %>%
    rename(obs = PAdata) %>%
    ungroup() %>%
    dplyr::left_join(site_data %>%
                         dplyr::select(Name, site) %>%
                         distinct(),
                     by = c("SiteName" = "Name")) %>%
    left_join(site_data %>%
                  dplyr::select(year, occasion) %>%
                  distinct(),
              by = "year") %>%
    group_by(site, occasion) %>%
    mutate(visit = row_number()) %>%
    ungroup() %>%
    data.table::as.data.table()


# Smooth for spatial effect on psi
fit <- fit_occu(list(psi ~ 1 + prcp + log(water+0.1) +
                         t2(lon, lat, occasion, k = c(25, 3), bs = c("ts", "cs"), d = c(2,1)),
                     p ~ 1 + log(TotalHours+0.1) + s(month, bs = "cs")),
                visit_data, site_data)

saveRDS(fit, paste0("analysis/output/", sp_sel, "_occur_fit.rds"))


# Explore model -----------------------------------------------------------

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

pred <- predict(fit, visit_data,  pred_data, nboot = 1000)

pred_data %>%
    as.data.frame() %>%
    left_join(gm, by = "Name") %>%
    sf::st_sf() %>%
    mutate(psi = pred$psi) %>%
    ggplot() +
    geom_sf(aes(fill = psi), size = 0.01) +
    scale_fill_viridis_c() +
    facet_wrap("year")

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
det <- visit_data %>%
    group_by(SiteName, year) %>%
    summarise(dets = sum(obs)) %>%
    ungroup() %>%
    mutate(det = if_else(dets > 0, 1L, 0L))

gm %>%
    right_join(det, by = c("Name" = "SiteName")) %>%
    ggplot() +
    geom_sf(aes(fill = factor(det)), size = 0.01) +
    scale_fill_viridis_d(option = "E") +
    facet_wrap("year")


# Explore non-linear effects ----------------------------------------------

pred_data <- pred_data %>%
    mutate(prcp = min(prcp),
           water = min(water)) %>%
    data.table::as.data.table()

pred <- predict(fit, visit_data,  pred_data, nboot = 0)

pred_data %>%
    as.data.frame() %>%
    left_join(gm, by = "Name") %>%
    sf::st_sf() %>%
    mutate(psi = pred$psi) %>%
    ggplot() +
    geom_sf(aes(fill = psi), size = 0.01) +
    scale_fill_viridis_c() +
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
