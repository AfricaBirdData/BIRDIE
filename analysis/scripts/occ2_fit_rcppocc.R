# 25-05-2021

# In this script we fit a spatial occupancy model to SABAP2 data

library(dplyr)
library(Rcppocc)
library(BIRDIE)

rm(list = ls())


# Select a species and a region -------------------------------------------

sp_sel <- 6

region <- "North West"


# Load occupancy data -----------------------------------------------------

occudata <- readRDS(paste0("analysis/data/pa_dat_", sp_sel, "_wcovts.rds"))


# Prepare spatial occupancy data ------------------------------------------

# Save pentads for later
ptd <- occudata$sites

# Subset one year
yr <- 2008

occudata$visits <- occudata$visits %>%
    filter(year == yr) %>%
    as.data.frame()

occudata$sites <- occudata$sites %>%
    sf::st_drop_geometry() %>%
    dplyr::select(Name, lon, lat, contains(as.character(yr)), water) %>%
    rename_with(.fn = ~gsub(paste0("_", yr), "", .x)) %>%
    as.data.frame()

# det is the difference of actual and potential evapotranspiration (aet - pet)
occudata$visits <- occudata$visits %>%
    dplyr::mutate(det = aet - pet,
                  tdif = tmax - tmin)

occudata$sites <- occudata$sites %>%
    dplyr::mutate(det = aet - pet,
                  tdif = tmax - tmin)

# Take log of water occurrence
occudata$sites <- occudata$sites %>%
    dplyr::mutate(det = aet - pet,
                  tdif = tmax - tmin,
                  log_water = log(water + 0.1))

# Define visit and site covariates
covts_s <- c("prcp", "tmax", "tdif", "log_water")
covts_v <- c("prcp", "tmax", "tdif")

# Standardize covariates
occudata$visits <- occudata$visits %>%
    dplyr::mutate(across(all_of(covts_v), ~scale(.x)))

occudata$sites <- occudata$sites %>%
    dplyr::mutate(across(all_of(covts_s), ~scale(.x)))

# Prepare data for stocc
so_data <- stocc::make.so.data(visit.data = occudata$visits,
                               site.data = occudata$sites,
                               names = list(visit = list(site = "SiteName", obs = "PAdata"),
                                            site = list(site = "Name", coords = c("lon", "lat"))))


# Fit simple model --------------------------------------------------------

# Set the prior distributions used
nalphas <- 2
nbetas <- length(covts_s)+1

beta_m <- matrix(rep(0, nbetas), ncol = 1)
sigma_inv_beta_p <- diag(nbetas)/1 # prior inverse covariance for beta
alpha_m <- matrix(rep(0, nalphas), ncol = 1)
sigma_inv_alpha_p <- diag(nalphas)/1 # prior inverse covariance for alpha

st <- Sys.time()
fit <- occSPATlogit(detection.model = reformulate(c("log(TotalHours)")),
                    occupancy.model = reformulate(covts_s, "log(water+0.1)"),
                    spatial.model = list(threshold = 1,
                                         moran.cut = 0.1*nrow(so_data$site)),
                    so.data = so_data,
                    prior = list(a.tau = 0.5, b.tau = 0.0005, tau = 1,
                                 Q.d=sigma_inv_alpha_p, mu.d = alpha_m,
                                 Q.o=sigma_inv_beta_p, mu.o = beta_m),
                    control = list(ndraws = 1000, percent_burn_in = 0.5))
en <- Sys.time()
en - st

future::plan("multisession", workers = 3)
fit <- BIRDIE::prllOccSPATlogit(nchains = 3,
                                detection.model = reformulate(c("log(TotalHours)")),
                                occupancy.model = reformulate(covts_s, "log(water+0.1)"),
                                spatial.model = list(threshold = 1,
                                                     moran.cut = 0.1*nrow(so_data$site)),
                                so.data = so_data,
                                prior = list(a.tau = 0.5, b.tau = 0.0005, tau = 1,
                                             Q.d=sigma_inv_alpha_p, mu.d = alpha_m,
                                             Q.o=sigma_inv_beta_p, mu.o = beta_m),
                                control = list(ndraws = 1000, percent_burn_in = 0.5))
future::plan("sequential")


saveRDS(fit, paste0("analysis/output/", sp_sel, "_so_fit.rds"))


# Explore model multichain -----------------------------------------------------

library(ggplot2)

fit <- readRDS(paste0("analysis/output/", sp_sel, "_so_fit.rds"))

# Plot detection parameters
extractParRcppocc(fit, "alpha") %>%
    setNames(c("Intcp", "log_hours", "chain")) %>%
    group_by(chain) %>%
    mutate(iter = row_number()) %>%
    tidyr::pivot_longer(cols = -c(iter, chain), names_to = "param", values_to = "values") %>%
    ggplot() +
    geom_line(aes(x = iter, y = values, group = chain, col = factor(chain))) +
    facet_grid(param~.)


# Plot occupancy parameters
extractParRcppocc(fit, "beta") %>%
    setNames(c("Intercept", covts_s, "chain")) %>%
    group_by(chain) %>%
    mutate(iter = row_number()) %>%
    tidyr::pivot_longer(cols = -c(iter, chain), names_to = "param", values_to = "values") %>%
    ggplot() +
    geom_line(aes(x = iter, y = values, group = chain, col = factor(chain))) +
    facet_grid(param~.)


# Occupancy probabilities
ptd %>%
    mutate(real.occ = extractParRcppocc(fit, "real.occ") %>%
               dplyr::select(-chain) %>%
               apply(., 2, mean)) %>%
    ggplot() +
    geom_sf(aes(fill = real.occ)) +
    scale_fill_viridis_c(limits = c(0,1))

ptd %>%
    mutate(real.occ = extractParRcppocc(fit, "psi") %>%
               dplyr::select(-chain) %>%
               apply(., 2, mean)) %>%
    ggplot() +
    geom_sf(aes(fill = real.occ)) +
    scale_fill_viridis_c(limits = c(0,1))

# Plot covariates
ptd %>%
    dplyr::select(Name, lon, lat, contains(as.character(yr)), water) %>%
    ggplot() +
    geom_sf(aes(fill = log(water))) +
    scale_fill_viridis_c()

ptd %>%
    dplyr::select(Name, lon, lat, contains(as.character(yr)), water) %>%
    ggplot() +
    geom_sf(aes(fill = prcp_2008)) +
    scale_fill_viridis_c()

ptd %>%
    dplyr::select(Name, lon, lat, contains(as.character(yr)), water) %>%
    ggplot() +
    geom_sf(aes(fill = tmax_2008)) +
    scale_fill_viridis_c()

ptd %>%
    dplyr::select(Name, lon, lat, contains(as.character(yr)), water) %>%
    ggplot() +
    geom_sf(aes(fill = tmax_2008 - tmin_2008)) +
    scale_fill_viridis_c()


# Explore model one chain ------------------------------------------------------

library(ggplot2)

fit <- readRDS(paste0("analysis/output/", sp_sel, "_so_fit.rds"))

# Plot detection parameters
# fit$alpha %>%
#     t() %>%
# as.data.frame() %>%

extractParRcppocc(fit, "alpha") %>%
    setNames(c("Intcp", "log_hours", "chain")) %>%
    group_by(chain) %>%
    mutate(iter = row_number()) %>%
    tidyr::pivot_longer(cols = -c(iter, chain), names_to = "param", values_to = "values") %>%
    ggplot() +
    geom_line(aes(x = iter, y = values, group = chain, col = factor(chain))) +
    facet_grid(param~.)


# Plot occupancy parameters
# fit$beta %>%
#     t() %>%
#     as.data.frame() %>%

extractParRcppocc(fit, "beta") %>%
    setNames(c("Intercept", covts_s, "chain")) %>%
    group_by(chain) %>%
    mutate(iter = row_number()) %>%
    tidyr::pivot_longer(cols = -c(iter, chain), names_to = "param", values_to = "values") %>%
    ggplot() +
    geom_line(aes(x = iter, y = values, group = chain, col = factor(chain))) +
    facet_grid(param~.)


# Occupancy probabilities
ptd %>%
    mutate(real.occ = apply(fit$real.occ, 1, mean)) %>%
    ggplot() +
    geom_sf(aes(fill = real.occ)) +
    scale_fill_viridis_c(limits = c(0,1))

ptd %>%
    mutate(psi = apply(fit$psi, 1, mean)) %>%
    ggplot() +
    geom_sf(aes(fill = psi)) +
    scale_fill_viridis_c(limits = c(0,1))

# Plot covariates
ptd %>%
    dplyr::select(Name, lon, lat, contains(as.character(yr)), water) %>%
    ggplot() +
    geom_sf(aes(fill = log(water))) +
    scale_fill_viridis_c()
