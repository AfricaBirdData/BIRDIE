# 25-05-2021

# In this script we fit a spatial occupancy model to SABAP2 data

library(tidyverse)
library(stocc)
library(BIRDIE)

rm(list = ls())


# Load occupancy data -----------------------------------------------------

occudata <- readRDS("analysis/data/occudata.rds")

# Load covariates
covts <- readRDS("analysis/data/pentd_nw_prcp.rds")


# Prepare spatial occupancy data ------------------------------------------

# Subset one year
yr <- 2008

# Prepare visit data
visit_data <- occudata %>%
    dplyr::filter(year %in% yr) %>%
    dplyr::select(Pentad, lon, lat, year, TotalHours, detc, StartDate) %>%
    as.data.frame()

# Add basis values of non-linear effect for time of year
visit_data <- visit_data %>%
    mutate(month = lubridate::month(StartDate))

spl_bs <- cSplineDes(x = visit_data$month,
                     knots = 1:12,
                     ord = 4, derivs = 0)
colnames(spl_bs) <- paste0("cyclic.", 1:ncol(spl_bs))
visit_data <- cbind(visit_data, spl_bs)

# Prepare site data
site_data <- covts %>%
    sf::st_drop_geometry() %>%
    dplyr::select(Pentad = Name, lon, lat, avg_prcp = contains(as.character(yr))) %>%
    mutate(intcp = 1)
# %>%
#     pivot_longer(cols = contains("avg_prcp"), names_to = "year", values_to = "avg_prcp") %>%
#     as.data.frame()

so_data <- make.so.data(visit.data = visit_data,
                        site.data = site_data,
                        names = list(visit = list(site = "Pentad", obs = "detc"),
                                     site = list(site = "Pentad", coords = c("lon", "lat"))))


# Fit occupancy model -----------------------------------------------------
reformulate(c("TotalHours", paste0("cyclic.", 1:(ncol(spl_bs)))))
# This spatial model does not allow the inclusion of year as a covariate for occupancy
fit <- spatial.occupancy(
    detection.model = reformulate(c("log(TotalHours)", paste0("cyclic.", 1:(ncol(spl_bs))))),
    # detection.model = ~log(TotalHours),
    occupancy.model = ~ 1 + avg_prcp,
    spatial.model = list(model="rsr", threshold=1, moran.cut=0.1*nrow(site_data)),
    so.data = so_data,
    prior = list(a.tau=0.5, b.tau=0.00005,
                 mu.b = 0.0, Q.b=0.1,
                 mu.g = 0.0, Q.g=0.1),
    control = list(burnin=1000/5, iter=4000/5, thin=5)
)

# Save fit
saveRDS(fit, "analysis/output/so_fit.rds")



# Explore model -----------------------------------------------------------

fit <- readRDS("analysis/output/so_fit.rds")

# Plot detection parameters
plot(fit$beta)

# Plot occupancy parameters
plot(fit$gamma)

# Occupancy probabilities
fitdf <- fit$occupancy.df

fitdf %>%
    ggplot() +
    geom_point(aes(x = lon, y = lat, colour = psi.est)) +
    scale_colour_viridis_c() +
    coord_equal()

fit$occupancy.df %>%
    ggplot() +
    geom_point(aes(x = lon, y = lat, col = eta.est)) +
    scale_colour_viridis_c() +
    coord_equal()

fit$occupancy.df %>%
    ggplot() +
    geom_point(aes(x = lon, y = lat, col = real.occ.est)) +
    scale_colour_viridis_c() +
    coord_equal()
