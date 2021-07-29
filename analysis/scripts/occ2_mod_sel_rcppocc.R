# 25-05-2021

# In this script we fit a spatial occupancy model to SABAP2 data

library(dplyr)
library(Rcppocc)
library(BIRDIE)
library(ggplot2)
library(tidyr)

rm(list = ls())


# Select a species and a region -------------------------------------------

sp_sel <- 41

region <- "North West"


# Load occupancy data -----------------------------------------------------

occudata <- readRDS(paste0("analysis/data/pa_dat_", sp_sel, "_wcovts_nw.rds"))

# Save pentads spatial object for later
pentads <- occudata$sites %>%
    dplyr::select(Name, geometry)

# Prepare spatial occupancy data ------------------------------------------

# Subset one year
yr <- 2008

occudata$visits <- occudata$visits %>%
    filter(year == yr) %>%
    as.data.frame()

occudata$sites <- occudata$sites %>%
    sf::st_drop_geometry() %>%
    dplyr::select(Name, lon, lat, contains(as.character(yr))) %>%
    rename_with(.fn = ~gsub(paste0("_", yr), "", .x))

# det is the difference of actual and potential evapotranspiration (aet - pet)
occudata$visits <- occudata$visits %>%
    dplyr::mutate(det = aet - pet,
                  tdif = tmax - tmin,
                  log_hours = log(TotalHours))

occudata$sites <- occudata$sites %>%
    dplyr::mutate(det = aet - pet,
                  tdif = tmax - tmin)

# Define covariate names
covts <- c("prcp", "tmax", "tmin", "tdif", "aet", "pet", "det")

# Standardize covariates
occudata$visits <- occudata$visits %>%
    dplyr::mutate(across(all_of(covts), ~scale(.x)))

occudata$sites <- occudata$sites %>%
    dplyr::mutate(across(all_of(covts), ~scale(.x)))

# Prepare data for stocc
so_data <- stocc::make.so.data(visit.data = occudata$visits,
                               site.data = occudata$sites,
                               names = list(visit = list(site = "SiteName", obs = "PAdata"),
                                            site = list(site = "Name", coords = c("lon", "lat"))))


# Fit model ---------------------------------------------------------------

# Define visit and site covariates
covts_v <- c("log_hours", paste0("toy.", 1:4))
covts_s <- c("prcp", "tmax", "tdif")

# Set the prior distributions used
nalphas <- length(covts_v)
nbetas <- length(covts_s) + 1

beta_m <- matrix(rep(0, nbetas), ncol = 1)
sigma_inv_beta_p <- diag(nbetas)/1 # prior inverse covariance for beta
alpha_m <- matrix(rep(0, nalphas), ncol = 1)
sigma_inv_alpha_p <- diag(nalphas)/1 # prior inverse covariance for alpha

fit <- occSPATlogit(detection.model = reformulate(c(-1, covts_v)),
                    occupancy.model = reformulate(covts_s),
                    spatial.model = list(threshold = 1,
                                         moran.cut = 0.1*nrow(so_data$site)),
                    so.data = so_data,
                    prior = list(a.tau = 0.5, b.tau = 0.0005, tau = 1,
                                 Q.d=sigma_inv_alpha_p, mu.d = alpha_m,
                                 Q.o=sigma_inv_beta_p, mu.o = beta_m),
                    control = list(ndraws = 4000, percent_burn_in = 0.5))

# Save fit
saveRDS(fit, paste0("analysis/output/", sp_sel, "_so_fit.rds"))


# Explore model -----------------------------------------------------------

fit <- readRDS(paste0("analysis/output/", sp_sel, "_so_fit.rds"))

# Plot detection parameters
fit$alpha %>%
    t() %>%
    as.data.frame() %>%
    setNames(covts_v) %>%
    mutate(iter = row_number()) %>%
    pivot_longer(cols = -iter, names_to = "param", values_to = "values") %>%
    ggplot() +
    geom_line(aes(x = iter, y = values)) +
    facet_grid(param~.)


# Plot occupancy parameters
fit$beta %>%
    t() %>%
    as.data.frame() %>%
    setNames(c("Intercept", covts_s)) %>%
    mutate(iter = row_number()) %>%
    pivot_longer(cols = -iter, names_to = "param", values_to = "values") %>%
    ggplot() +
    geom_line(aes(x = iter, y = values)) +
    facet_grid(param~.)

# Occupancy probabilities
pentads %>%
    cbind(data.frame(psi.est = rowMeans(fit$psi),
           psi.sd = apply(fit$psi, 1, sd))) %>%
    ggplot() +
    geom_sf(aes(fill = psi.est), lwd = 0.1) +
    scale_fill_viridis_c(limits = c(0,1))

# Estimated occupancy
pentads %>%
    cbind(data.frame(occ.est = rowMeans(fit$real.occ),
                     occ.sd = apply(fit$real.occ, 1, sd))) %>%
    ggplot() +
    geom_sf(aes(fill = occ.est), lwd = 0.1) +
    scale_fill_viridis_c(limits = c(0,1))

# Estimated vs Observed occupancy
pentads_det <- occudata$visits %>%
    group_by(SiteName) %>%
    summarise(det = sum(PAdata)) %>%
    mutate(det = if_else(det > 0, 1L, 0L))

pentads %>%
    cbind(data.frame(occ.est = rowMeans(fit$real.occ),
                     occ.sd = apply(fit$real.occ, 1, sd))) %>%
    left_join(pentads_det[,c("SiteName", "det")], by = c("Name" = "SiteName")) %>%
    ggplot() +
    geom_sf(aes(fill = occ.est), size = 0.1) +
    geom_sf(data = . %>% filter(!is.na(det)), aes(colour = factor(det)),
            fill = NA, size = 0.3) +
    scale_fill_viridis_c(limits = c(0,1)) +
    scale_colour_viridis_d(option = "B")

pentads %>%
    cbind(data.frame(occ.est = rowMeans(fit$real.occ),
                     occ.sd = apply(fit$real.occ, 1, sd))) %>%
    left_join(pentads_det[,c("SiteName", "det")], by = c("Name" = "SiteName")) %>%
    ggplot() +
    geom_sf(aes(fill = occ.sd), size = 0.1) +
    geom_sf(data = . %>% filter(!is.na(det)), aes(colour = factor(det)),
            fill = NA, size = 0.3) +
    scale_fill_viridis_c(limits = c(0,1)) +
    scale_colour_viridis_d(option = "B")
