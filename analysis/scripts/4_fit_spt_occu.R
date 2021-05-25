# 25-05-2021

# In this script we fit a spatial occupancy model to SABAP2 data

library(tidyverse)
library(stocc)
devtools::load_all() # To load BIRDIE package

rm(list = ls())


# Load data ---------------------------------------------------------------

occudata <- readRDS("analysis/data/occudata.rds")


# Prepare spatial occupancy data ------------------------------------------

# Prepare visit data
visit_data <- occudata %>%
    dplyr::select(Pentad, lon, lat, year, TotalHours, detc) %>%
    as.data.frame()

# Prepare site data
site_data <- occudata %>%
    group_by(Pentad) %>%
    summarize(lon = first(lon),
              lat = first(lat),
              intcp = 1) %>%
    as.data.frame()

so_data <- make.so.data(visit.data = visit_data,
                        site.data = site_data,
                        names = list(visit = list(site = "Pentad", obs = "detc"),
                                     site = list(site = "Pentad", coords = c("lon", "lat"))))


# Fit occupancy model -----------------------------------------------------

# This spatial model does not allow the inclusion of year as a covariate for occupancy
fit <- spatial.occupancy(
    detection.model = ~ TotalHours,
    occupancy.model = ~ 1,
    spatial.model = list(model = "icar", threshold = 1),
    so.data = so_data,
    prior = list(a.tau=0.5, b.tau=0.00005, Q.b=0.1, Q.g=0.1),
    control = list(burnin=1000/5, iter=4000/5, thin=5)
)

# Save fit
saveRDS(fit, "analysis/output/so_fit.rds")

# Plot detection parameters
plot(fit$beta)

# Plot occupancy parameters
plot(fit$gamma)

# Occupancy probabilities
fit$occupancy.df

fit$occupancy.df %>%
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = psi.est)) +
    scale_fill_viridis_c() +
    coord_equal()

fit$occupancy.df %>%
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = eta.est)) +
    scale_fill_viridis_c() +
    coord_equal()

fit$occupancy.df %>%
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = real.occ.est)) +
    scale_fill_viridis_c() +
    coord_equal()
