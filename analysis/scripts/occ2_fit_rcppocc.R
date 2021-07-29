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

occudata <- readRDS(paste0("analysis/data/pa_dat_", sp_sel, "_wcovts_nw.rds"))


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

fit <- occSPATlogit(detection.model = reformulate(c("log(TotalHours)")),
                    occupancy.model = reformulate(covts_s, "log(water+0.1)"),
                    spatial.model = list(threshold = 1,
                                         moran.cut = 0.1*nrow(so_data$site)),
                    so.data = so_data,
                    prior = list(a.tau = 0.5, b.tau = 0.0005, tau = 1,
                                 Q.d=sigma_inv_alpha_p, mu.d = alpha_m,
                                 Q.o=sigma_inv_beta_p, mu.o = beta_m),
                    control = list(ndraws = 1000, percent_burn_in = 0.5))

saveRDS(fit, paste0("analysis/output/", sp_sel, "_so_fit.rds"))

# Explore model -----------------------------------------------------------

fit <- readRDS(paste0("analysis/output/", sp_sel, "_so_fit.rds"))

# Plot detection parameters
fit$alpha %>%
    t() %>%
    as.data.frame() %>%
    setNames(c("Intcp", "log_hours")) %>%
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
library(ggplot2)

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

fit %>%
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

occuobs <- occudata %>%
    group_by(Pentad, detc) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    tidyr::complete(Pentad, tidyr::nesting(detc), fill = list(n = 0)) %>%
    tidyr::pivot_wider(names_from = detc, values_from = n, names_prefix = "det") %>%
    mutate(effort = det0 + det1,
           reprate = det1/effort) %>%
    left_join(dplyr::select(occudata, Pentad, lon, lat), by = "Pentad")

occuobs %>%
    ggplot() +
    geom_point(aes(x = lon, y = lat, col = reprate)) +
    scale_colour_viridis_c() +
    coord_equal()


library(Rcppocc)

# Set the prior distributions used
nalphas <- 13
nbetas <- 2

beta_m<-matrix( rep(0,nbetas), ncol=1)
sigma_inv_beta_p<-diag(nbetas)/1000 #prior inverse covariance for beta
alpha_m<-matrix( rep(0,nalphas), ncol=1)
sigma_inv_alpha_p<-diag(nalphas)/1000 #prior inverse covariance for alpha

# Rename properly
so_data$visit$PAdata <- so_data$visit$detc
so_data$visit$SiteName <- so_data$visit$Pentad

fit <- occSPATlogit(detection.model = reformulate(c("1 + log(TotalHours)", paste0("cyclic.", 1:11))),
                    # detection.model = ~log(TotalHours),
                    occupancy.model = reformulate(c(1, paste("prcp", yr, sep = "_"))),
                    spatial.model = list(threshold = 1,
                                         moran.cut = 0.1*nrow(so_data$site)),
                    so.data = so_data,
                    prior = list(a.tau = 0.5, b.tau = 0.0005,
                                 tau=1,
                                 Q.d=sigma_inv_alpha_p, mu.d = alpha_m,
                                 Q.o=sigma_inv_beta_p, mu.o = beta_m),
                    control = list(ndraws = 1000, percent_burn_in = 0.5))

apply(fit$real.occ, 1, mean)


# Model selection ---------------------------------------------------------

mods <- list("log(TotalHours)",
             c("log(TotalHours)", covts_v),
             c("log(TotalHours)", paste0("toy.", 1:4)),
             c(covts_v, "log(TotalHours)", paste0("toy.", 1:4)))

PPL <- vector("list", length = length(mods))

for(i in seq_along(mods)){

    # Set the prior distributions used
    nalphas <- length(mods[[i]])
    nbetas <- length(covts_s)

    beta_m <- matrix(rep(0, nbetas), ncol = 1)
    sigma_inv_beta_p <- diag(nbetas)/1000 # prior inverse covariance for beta
    alpha_m <- matrix(rep(0, nalphas), ncol = 1)
    sigma_inv_alpha_p <- diag(nalphas)/1000 # prior inverse covariance for alpha

    fit <- occSPATlogit(detection.model = reformulate(mods[[i]]),
                        occupancy.model = reformulate(covts_s),
                        spatial.model = list(threshold = 1,
                                             moran.cut = 0.1*nrow(so_data$site)),
                        so.data = so_data,
                        prior = list(a.tau = 0.5, b.tau = 0.0005, tau = 1,
                                     Q.d=sigma_inv_alpha_p, mu.d = alpha_m,
                                     Q.o=sigma_inv_beta_p, mu.o = beta_m),
                        control = list(ndraws = 1000, percent_burn_in = 0.5))

    PPL[[i]] <- unlist(fit[c("D.m", "G.m", "P.m")])

}

PPL


# Fit best model ----------------------------------------------------------

# Set the prior distributions used
nalphas <- length(mods[[1]])
nbetas <- length(covts_s) + 1

beta_m <- matrix(rep(0, nbetas), ncol = 1)
sigma_inv_beta_p <- diag(nbetas)/1 # prior inverse covariance for beta
alpha_m <- matrix(rep(0, nalphas), ncol = 1)
sigma_inv_alpha_p <- diag(nalphas)/1 # prior inverse covariance for alpha

fit <- occSPATlogit(detection.model = reformulate(c(-1, mods[[1]])),
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

