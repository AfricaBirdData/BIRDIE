library(occuR)
library(tidyverse)
library(BIRDIE)

rm(list = ls())


# Load general data -------------------------------------------------------

# We want to select species present at Barberspan
bbpan <- barberspan

# Sample 10 species
set.seed(378465)
spp <- sample(unique(bbpan$spp), 10)
set.seed(NULL)

# Load site data
sitedata <- readRDS("analysis/data/site_dat_sa_wcovts.rds")


# Define models -----------------------------------------------------------

# Detection.
# I wanted to test whether a smooth effect of TotalHours works better than
# the log of TotalHours. However, we run into fitting problems when TotalHours
# are very high (> 100 or so).
visit_covts <- list(mod1 = "1",
                    mod2 = c("1", "log(TotalHours+1)"),
                    # mod3 = c("1", "s(TotalHours, bs = 'cs')"),
                    mod3 = c("1", "log(TotalHours+1)", "s(month, bs = 'cs')"))

# Occupancy
site_covts <- list(mod1 = "1",
                   mod2 = c("1", "log(water+0.1)"),
                   mod3 = c("1", "s(water, bs = 'cs')"),
                   mod4 = c("1", "s(water, bs = 'cs')", "prcp"))

# Create output lists
aic_site <- vector("list", length = length(site_covts))
aic_visit <- vector("list", length = length(visit_covts))


# Iterate and fit models --------------------------------------------------

for(i in seq_along(spp)){

    # Prepare visit data
    sp_sel <- spp[i]

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


    # Fit detection models -----------------------------------------------------

    future::plan("multisession", workers = 6)

    aic_visit[[i]] <- furrr::future_map2_dfr(visit_covts, names(visit_covts),
                                             ~selOccuRmod(forms = list(reformulate(.x, response = "p"),
                                                                       reformulate(site_covts[[1]], response = "psi")),
                                                          type = "visit",
                                                          visit_data = occuRdata$visit,
                                                          site_data = occuRdata$site,
                                                          mod_id = paste0("sp_", sp_sel, "_", .y)),
                                             .options = furrr::furrr_options(packages = "occuR"))


    # Select site model ------------------------------------------------------

    aic_site[[i]] <- furrr::future_map2_dfr(site_covts, names(site_covts),
                                            ~selOccuRmod(forms = list(reformulate(visit_covts[[1]], response = "p"),
                                                                      reformulate(.x, response = "psi")),
                                                         type = "site",
                                                         visit_data = occuRdata$visit,
                                                         site_data = occuRdata$site,
                                                         mod_id = paste0("sp_", sp_sel, "_", .y)),
                                            .options = furrr::furrr_options(packages = "occuR"))

    future::plan("sequential")

}



# Create dataframes
aic_visit <- do.call(rbind, aic_visit)
aic_site <- do.call(rbind, aic_site)

# Save results
saveRDS(aic_visit, file = "analysis/output/aic_visit.rds")
saveRDS(aic_site, file = "analysis/output/aic_site.rds")

# Calculate Akaike weights
aic_visit <- aic_visit %>%
    group_by(species) %>%
    mutate(delta_aic = AIC - min(AIC),
           lik_aic = exp(-0.5*delta_aic),
           w = round(lik_aic/sum(lik_aic), 3)) %>%
    ungroup()

aic_site <- aic_site %>%
    group_by(species) %>%
    mutate(delta_aic = AIC - min(AIC),
           lik_aic = exp(-0.5*delta_aic),
           w = round(lik_aic/sum(lik_aic), 3)) %>%
    ungroup()

# Check that weights add up to one
aic_visit %>%
    group_by(species) %>%
    summarize(total = sum(w))




# Plot
aic_site %>%
    ggplot() +
    geom_boxplot(aes(x = factor(mod), y = w)) +
    geom_jitter(aes(x = factor(mod), y = w), col = "red", alpha = 0.5)

aic_visit %>%
    ggplot() +
    geom_boxplot(aes(x = factor(mod), y = w)) +
    geom_jitter(aes(x = factor(mod), y = w), col = "red", alpha = 0.5)



# Test new models ---------------------------------------------------------

# Load existing AIC scores
aic_visit <- readRDS(file = "analysis/output/aic_visit.rds")
aic_site <- readRDS(file = "analysis/output/aic_site.rds")


# Iterate and fit models --------------------------------------------------

for(i in seq_along(spp)){

    # Prepare visit data
    sp_sel <- spp[i]

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

    # Define new models
    new_visit_mod <- list(c("1", "log(TotalHours+1)", "s(prcp, bs = 'cs')"))
    new_site_mod <- list(c("1", "s(water, bs = 'cs')", "s(prcp, bs = 'cs')"))

    # Maximum model index
    max_visit_idx <- max(as.numeric(gsub("mod", "", unique(aic_visit$mod))))
    max_site_idx <- max(as.numeric(gsub("mod", "", unique(aic_site$mod))))

    # Name new models
    names(new_visit_mod) <- paste0("mod", max_visit_idx+1)
    names(new_site_mod) <- paste0("mod", max_site_idx+1)

    # Fit new visit model
    fit <- fit_occu(forms = list(reformulate(new_visit_mod[[1]], response = "p"),
                                 reformulate("1", response = "psi")),
                    visit_data = occuRdata$visit,
                    site_data = occuRdata$site)

    new_aic_visit <- aics %>%
        mutate(form = map(visit_covts, ~reformulate(.x)),
               species = sp_sel,
               mod = names(visit_covts))

    # Fit new site model
    fit <- fit_occu(forms = list(reformulate("1", response = "p"),
                             reformulate(new_site_mod[[1]], response = "psi")),
                    visit_data = occuRdata$visit,
                    site_data = occuRdata$site)

    new_aic_site <- data.frame(df = dof.occuR(fit, each = FALSE),
                               form = new_site_mod[[1]],
                               AIC = AIC(fit),
                               species = sp_sel[i],
                               mod = names(new_site_mod),
                               delta_aic = NA,
                               lik_aic = NA,
                               w = NA)
