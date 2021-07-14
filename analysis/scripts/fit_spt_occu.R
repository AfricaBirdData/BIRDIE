# 25-05-2021

# In this script we fit a spatial occupancy model to SABAP2 data

library(dplyr)
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

# Rename covariates
covts <- covts %>%
    rename_with(.fn = ~gsub("X", "prcp_", .x), .cols = contains("X"))

# Prepare data
so_data <- prepOccuData(occudata, yr, covts)


# Fit occupancy model -----------------------------------------------------

# This spatial model does not allow the inclusion of year as a covariate for occupancy
fit <- spatial.occupancy(
    detection.model = reformulate(c("log(TotalHours)", paste0("cyclic.", 1:11))),
    # detection.model = ~log(TotalHours),
    occupancy.model = reformulate(c(1, paste("prcp", yr, sep = "_"))),
    spatial.model = list(model="rsr", threshold = 1,
                         moran.cut = 0.1*nrow(so_data$site)),
    so.data = so_data,
    prior = list(a.tau=0.5, b.tau=0.00005,
                 Q.b=0.1,
                 Q.g=0.1),
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

library(ggplot2)

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



spat_logit <- occSPATlogit(detection.model = reformulate(c("1 + log(TotalHours)", paste0("cyclic.", 1:11))),
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

apply(spat_logit$real.occ, 1, mean)
