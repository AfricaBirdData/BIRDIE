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



# Define models -----------------------------------------------------------

visit_covts <- list(mod1 = "1",
                    mod2 = c("1", "log(TotalHours+1)"),
                    mod3 = c("1", "s(TotalHours, bs = 'cs')"),
                    mod4 = c("1", "s(TotalHours, bs = 'cs')", "s(month, bs = 'cs')"))

site_covts <- list(mod1 = "1",
                   mod2 = c("1", "log(water+0.1)"),
                   mod3 = c("1", "s(water, bs = 'cs')"),
                   mod4 = c("1", "s(water, bs = 'cs')", "prcp"))


# Select visit model ------------------------------------------------------

fits_visit <- vector("list", length = length(visit_covts))

for(i in seq_along(visit_covts)){

    fit <- fit_occu(forms = list(reformulate(visit_covts[[i]], response = "p"),
                                 reformulate(site_covts[[1]], response = "psi")),
                    visit_data = visit_data,
                    site_data = site_data)

    fits_visit[[i]] <- fit

}

#do.call("AIC", fits_visit) # kills the session for some reason

aic_visit <- AIC(fits_visit[[1]], fits_visit[[2]], fits_visit[[3]], fits_visit[[4]])

aic_visit <- aic_visit %>%
    mutate(form = map(visit_covts, ~reformulate(.x)))

saveRDS(aic_visit, file = "analysis/output/aic_visit.rds")


# Select site model ------------------------------------------------------

fits_site <- vector("list", length = length(site_covts))

for(i in seq_along(site_covts)){

    fit <- fit_occu(forms = list(reformulate(visit_covts[[1]], response = "p"),
                                 reformulate(site_covts[[i]], response = "psi")),
                    visit_data = visit_data,
                    site_data = site_data)

    fits_site[[i]] <- fit

}

#do.call("AIC", fits) # kills the session for some reason

aic_site <- AIC(fits_site[[1]], fits_site[[2]], fits_site[[3]], fits_site[[4]])

aic_site <- aic_site %>%
    mutate(form = map(site_covts, ~reformulate(.x)))

saveRDS(aic_site, file = "analysis/output/aic_site.rds")

