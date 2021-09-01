# 25-05-2021

# In this script we fit a spatial occupancy model to SABAP2 data

library(dplyr)
library(stocc)
library(BIRDIE)

rm(list = ls())


# Select a species and a region -------------------------------------------

sp_sel <- 41

region <- "North West"


# Load occupancy data -----------------------------------------------------

occudata <- readRDS(paste0("analysis/data/pa_dat_", sp_sel, "_wcovts_nw.rds"))


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

# Define visit and site covariates
covts_s <- c("prcp", "tmax")
covts_v <- c("prcp", "tmax")

# det is the difference of actual and potential evapotranspiration (aet - pet)
# Also, scale and centre covariates
# occudata$visits <- occudata$visits %>%
#     dplyr::mutate(det = aet - pet) %>%
#     dplyr::mutate(across(all_of(covts_v), ~scale(.x)))
#
# occudata$sites <- occudata$sites %>%
#     dplyr::mutate(det = aet - pet) %>%
#     dplyr::mutate(across(all_of(covts_s), ~scale(.x)))

# Prepare data for stocc
so_data <- stocc::make.so.data(visit.data = occudata$visits,
                               site.data = occudata$sites,
                               names = list(visit = list(site = "SiteName", obs = "PAdata"),
                                            site = list(site = "Name", coords = c("lon", "lat"))))


# Fit simple model --------------------------------------------------------

fit <- spatial.occupancy(
    detection.model = reformulate(c("log(TotalHours)")),
    occupancy.model = reformulate(covts_s),
    spatial.model = list(model="rsr", threshold = 1,
                         moran.cut = 0.1*nrow(so_data$site)),
    so.data = so_data,
    prior = list(a.tau=0.5, b.tau=0.00005,
                 Q.b=0.1,
                 Q.g=0.1),
    control = list(burnin=1000/5, iter=4000/5, thin=5)
)

unlist(fit[c("D.m", "G.m", "P.m")])


# Model selection ---------------------------------------------------------

mods <- list("log(TotalHours)",
             c("log(TotalHours)", covts_v),
             c("log(TotalHours)", paste0("toy.", 1:4)),
             c(covts_v, "log(TotalHours)", paste0("toy.", 1:4)))

PPL <- vector("list", length = length(mods))

for(i in seq_along(mods)){

    fit <- spatial.occupancy(
        detection.model = reformulate(mods[[i]]),
        occupancy.model = reformulate(covts_s),
        spatial.model = list(model="rsr", threshold = 1,
                             moran.cut = 0.1*nrow(so_data$site)),
        so.data = so_data,
        prior = list(a.tau=0.5, b.tau=0.00005,
                     Q.b=0.1,
                     Q.g=0.1),
        control = list(burnin=1000/5, iter=4000/5, thin=5)
    )

    PPL[[i]] <- unlist(fit[c("D.m", "G.m", "P.m")])

}

PPL


# Fit best model ----------------------------------------------------------

fit <- spatial.occupancy(
    detection.model = reformulate(mods[[3]]),
    occupancy.model = reformulate(covts_s),
    spatial.model = list(model="rsr", threshold = 1,
                         moran.cut = 0.1*nrow(so_data$site)),
    so.data = so_data,
    prior = list(a.tau=0.5, b.tau=0.00005,
                 Q.b=0.1,
                 Q.g=0.1),
    control = list(burnin=2000/5, iter=8000/5, thin=5)
)

unlist(fit[c("D.m", "G.m", "P.m")])

# Save fit
saveRDS(fit, paste0("analysis/output/", sp_sel, "_so_fit.rds"))


# Explore model -----------------------------------------------------------

fit <- readRDS(paste0("analysis/output/", sp_sel, "_so_fit.rds"))

# Plot detection parameters
plot(fit$beta)

coda::traceplot(fit$beta)

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
