library(occuR)
library(BIRDIE)
library(dplyr)
library(sf)
library(ggplot2)

rm(list = ls())

birdie_dir <- "analysis/"

# Load general data -------------------------------------------------------

# We want to select species present at Barberspan
bbpan <- barberspan

# Sample 10 species
set.seed(378465)
spp <- sample(unique(bbpan$spp), 30)
set.seed(NULL)

# Load site data and subset year
year_sel <- 2010
sitedata <- readRDS(paste0(birdie_dir, "data/site_dat_sa_gee_08_19.rds")) %>%
    dplyr::select(Name, lon, lat, watocc_ever, dist_coast, elev, ends_with(match = as.character(year_sel)))

# I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE
sitedata <- sitedata %>%
    tidyr::drop_na()


# Define models -----------------------------------------------------------

# Detection.
# I wanted to test whether a smooth effect of TotalHours works better than
# the log of TotalHours. However, we run into fitting problems when TotalHours
# are very high (> 100 or so).
visit_covts <- list(mod1 = "1",
                    mod2 = c("1", "log(TotalHours+1)"),
                    mod3 = c("1", "log(TotalHours+1)", "s(month, bs = 'cs')"))

# Occupancy
coast_covts <- c("1", "dist_coast")

land_mods <- list(mod1 = "1",
                  # mod2 = c("1", "s(dist_coast, bs = 'cs')"),    # This one beats 3,4,5,6
                  # mod3 = c("1", "s(watocc_ever, bs = 'cs')"),
                  # mod4 = c("1", "s(watext, bs = 'cs')"),
                  # mod5 = c("1", "s(watrec, bs = 'cs')"),
                  # mod6 = c("1", "t2(watrec, watext, bs = 'ts', d = 2)"),
                  # mod7 = c("1", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"), # This one beats 3,4,5,6
                  # mod8 = c("1", "dist_coast", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"),
                  # mod9 = c("1", "s(prcp, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"),  # this one beats 8
                  # mod10 = c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"), # this one beats 9
                  # mod11 = c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"), # this one beats 10
                  # mod12 = c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"), # this one beats 11
                  # mod13 = c("1", "dist_coast", "elev", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"), # this doesn't converge
                  mod14 = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec", "t2(lon, lat, bs = 'ts')"), # better than 12
                  mod15 = c("dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')", "t2(lon, lat, bs = 'ts')")) # this doesn't converge

coast_mods <- list(mod1 = "1",
                   mod2 = c("1", "dist_coast"),
                   mod3 = c("1", "watocc_ever")#,
                   # mod4 = c("1", "watext"),
                   # mod5 = c("1", "watrec"),
                   # mod6 = c("1", "watrec * watext"),
                   # mod7 = c("1", "watext", "watrec"),               # This one beats 3,4,5,6
                   # mod8 = c("1", "dist_coast", "watext", "watrec"), # this one beats 9
                   # mod9 = c("1", "prcp", "watext", "watrec"),
                   # mod10 = c("1", "dist_coast", "prcp", "watext", "watrec"), #this one beats 8
                   # mod11 = c("1", "dist_coast", "prcp", "tdiff", "watext", "watrec"),   # this one beats 10 and 12
                   # mod12 = c("1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec"),
                   # mod13 = c("1", "dist_coast", "elev", "ndvi", "watext", "watrec"),
                   # mod14 = c("dist_coast", "prcp", "tdiff", "watext", "watrec", "t2(lon, lat, bs = 'ts')") # as good as 13
                   )

# Create output lists
aic_site <- vector("list", length = length(spp))
aic_visit <- vector("list", length = length(spp))
coast <- vector(length = length(spp))


# Iterate and fit models --------------------------------------------------

for(i in seq_along(spp)){
# for(i in 17:length(spp)){
    # Prepare visit data
    sp_sel <- spp[i]

    # Download species detection
    print("Downloading from ABAP")
    sp_detect <- ABAP::getAbapData(.spp_code = sp_sel,
                                   .region_type = "country",
                                   .region = "South Africa",
                                   .years = year_sel)


    # Load occupancy data -----------------------------------------------------

    # Load visit data, subset years and add detections
    visitdata <- readRDS(paste0(birdie_dir, "data/visit_dat_sa_gee_08_19.rds")) %>%
        filter(year == year_sel) %>%
        dplyr::left_join(sp_detect %>%
                      dplyr::select(CardNo, obs = Spp) %>%
                      mutate(obs = if_else(obs == "-", 0, 1)),
                  by = "CardNo")


    # Format to occuR ---------------------------------------------------------

    occuRdata <- prepDataOccuR(site_data = sitedata %>%
                                   sf::st_drop_geometry() %>%
                                   gatherYearFromVars(vars = names(.)[-c(1:6)], sep = "_") %>%
                                   mutate(tdiff = tmmx - tmmn),
                               visit_data = visitdata,
                               scaling = list(visit = NULL,
                                              site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec", "elev")))

    # Remove data from missing pentads? (THIS SHOULDN'T HAPPEN WHEN PENTADS ARE DOWNLOADED FROM THE API)
    occuRdata$visit <- occuRdata$visit %>%
        dplyr::filter(!is.na(site))

    occuRdata$site <- occuRdata$site %>%
        filter(site %in% unique(occuRdata$visit$site))

    m1 <- fit_occu(forms = list(reformulate(visit_covts[[1]], response = "p"),
                                reformulate(coast_covts, response = "psi")),
                   visit_data = occuRdata$visit,
                   site_data = occuRdata$site,
                   print = TRUE)

    if(m1$fit$par[2] < -2){
        mods <- coast_mods
        coast[i] <- 1
    } else {
        mods <- land_mods
        coast[i] <- 0
    }

    # Select visit model -----------------------------------------------------

    future::plan("multisession", workers = min(6, length(mods)))

    aic_visit[[i]] <- furrr::future_map2_dfr(visit_covts, names(visit_covts),
                                             ~selOccuRmod(forms = list(reformulate(.x, response = "p"),
                                                                       reformulate(mods[[1]], response = "psi")),
                                                          type = "visit",
                                                          visit_data = occuRdata$visit,
                                                          site_data = occuRdata$site,
                                                          mod_id = paste0("sp_", sp_sel, "_", .y)),
                                             .options = furrr::furrr_options(packages = "occuR"))

    aic_visit[[i]]$coast <- coast[[i]]


    # Select site model ------------------------------------------------------

    aic_site[[i]] <- furrr::future_map2_dfr(mods, names(mods),
                                            ~selOccuRmod(forms = list(reformulate(visit_covts[[1]], response = "p"),
                                                                      reformulate(.x, response = "psi")),
                                                         type = "site",
                                                         visit_data = occuRdata$visit,
                                                         site_data = occuRdata$site,
                                                         mod_id = paste0("sp_", sp_sel, "_", .y)),
                                            .options = furrr::furrr_options(packages = "occuR"))

    # Test single models
    # fits <- vector("list", length = length(mods))
    # for(m in seq_along(fits)){
    #     fits[[m]] <- fit_occu(forms = list(reformulate(visit_covts[[1]], response = "p"),
    #                                        reformulate(mods[[m]], response = "psi")),
    #                           # reformulate(site_covts_lin[[m]], response = "psi")),
    #                           visit_data = occuRdata$visit,
    #                           site_data = occuRdata$site)
    #
    #     print(dof.occuR(fits[[m]]))
    # }

    # # Plot effects?
    # plotOccuVars(occuRdata, "dist_coast")
    # plotOccuVarEffect(fit = fits[[2]],
    #                   new_occu_data = list(site = data.table(dist_coast = seq(-3, 3, length.out = 100)),
    #                                        visit = data.table(occasion = rep(1, 100))),
    #                   var = "dist_coast",
    #                   nboot = 1000)
    #
    # # Plot detections
    # plotDetections(site_data = rename(sitedata, site = Name),
    #                visit_data = rename(visitdata, site = Pentad))


    aic_site[[i]]$coast <- coast[[i]]


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
    mutate(species = gsub("sp_|_mod.*", "", mod),
           mod = gsub("sp_.*_", "", mod)) %>%
    group_by(species, coast) %>%
    mutate(delta_aic = AIC - min(AIC),
           lik_aic = exp(-0.5*delta_aic),
           w = round(lik_aic/sum(lik_aic), 3)) %>%
    ungroup()

aic_site <- aic_site %>%
    mutate(species = gsub("sp_|_mod.*", "", mod),
           mod = gsub("sp_.*_", "", mod)) %>%
    group_by(species, coast) %>%
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
    geom_jitter(aes(x = factor(mod), y = w), col = "red", alpha = 0.5) +
    facet_wrap("coast")

aic_visit %>%
    ggplot() +
    geom_boxplot(aes(x = factor(mod), y = w)) +
    geom_jitter(aes(x = factor(mod), y = w), col = "red", alpha = 0.5) +
    facet_wrap("coast")



# TEST NEW MODELS ---------------------------------------------------------

# Load existing AIC scores
aic_visit <- readRDS(file = "analysis/output/aic_visit.rds")
aic_site <- readRDS(file = "analysis/output/aic_site.rds")


# Iterate and fit models --------------------------------------------------

# Define new models (REMEMBER TO CHANGE MODEL NAME IN THE LIST)
new_visit_mod <- list(mod4 = c("1", "log(TotalHours+1)", "s(aet, bs = 'cs')"))
# new_site_mod <- list(mod8 = c("1", "s(water, bs = 'cs')", "s(prcp, bs = 'cs')", "s(aet/pet, bs = 'cs')"))

for(i in seq_along(spp)){
# for(i in 7:length(spp)){

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

    # Fit new visit model
    new_aic_visit <- selOccuRmod(forms = list(reformulate(new_visit_mod[[1]], response = "p"),
                                              reformulate("1", response = "psi")),
                                 type = "visit",
                                 visit_data = occuRdata$visit,
                                 site_data = occuRdata$site,
                                 mod_id = paste0("sp_", sp_sel, "_", names(new_visit_mod)))

    aic_visit <- rbind(aic_visit,
                       new_aic_visit)

    # Fit new site model
    # new_aic_site <- selOccuRmod(forms = list(reformulate("1", response = "p"),
    #                                          reformulate(new_site_mod[[1]], response = "psi")),
    #                             type = "site",
    #                             visit_data = occuRdata$visit,
    #                             site_data = occuRdata$site,
    #                             mod_id = paste0("sp_", sp_sel, "_", names(new_site_mod)))
    #
    #
    # aic_site <- rbind(aic_site,
    #                   new_aic_site)

}

# Save results
saveRDS(aic_visit, file = "analysis/output/aic_visit.rds")
saveRDS(aic_site, file = "analysis/output/aic_site.rds")

aic_visit <- readRDS("analysis/output/aic_visit.rds")
aic_site <- readRDS("analysis/output/aic_site.rds")

# Calculate Akaike weights
aic_visit <- aic_visit %>%
    mutate(species = gsub("sp_|_mod.*", "", mod),
           mod = gsub("sp_.*_", "", mod)) %>%
    group_by(species) %>%
    mutate(delta_aic = AIC - min(AIC),
           lik_aic = exp(-0.5*delta_aic),
           w = round(lik_aic/sum(lik_aic), 3)) %>%
    ungroup()

aic_site <- aic_site %>%
    mutate(species = gsub("sp_|_mod.*", "", mod),
           mod = gsub("sp_.*_", "", mod)) %>%
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

