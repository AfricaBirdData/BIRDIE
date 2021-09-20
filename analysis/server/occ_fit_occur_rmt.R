# In this script we fit an occupancy model with occuR for all species

library(occuR)
library(dplyr)
library(BIRDIE)

rm(list = ls())

birdie_dir <- "/home/birdie/analysis/"

# For now, we want to select species present at Barberspan
bbpan <- BIRDIE::barberspan %>%
    pull(spp) %>%
    unique()

for(i in seq_along(bbpan)){

    # Select a species and a region -------------------------------------------

    sp_sel <- bbpan[i]

    sp_name <- BIRDIE::barberspan %>%
        dplyr::filter(spp == sp_sel) %>%
        mutate(name = paste(taxon.Common_species, taxon.Common_group)) %>%
        pull(name) %>%
        unique()


    # Load occupancy site data -------------------------------------------------

    # Load site data
    sitedata <- readRDS(paste0(birdie_dir, "data/site_dat_sa_wcovts_08_19.rds"))


    # Prepare occupancy visit data --------------------------------------------

    future::plan("multisession", workers = 6)

    visitdata <- prepOccVisitData(region = "South Africa",
                                  sites = sitedata,
                                  species = sp_sel,
                                  years = 2008:2016,
                                  clim_covts = c("prcp", "tmax", "tmin", "aet", "pet"),
                                  covts_dir = paste0(birdie_dir, "downloads/"),
                                  file_fix = c("terraClim_", "_03_19"),
                                  savedir = paste0(birdie_dir, "data/pentads_sa.rds"))

    future::plan("sequential")


    saveRDS(visitdata, paste0(birdie_dir, "data/visit_dat_", sp_sel, "_wcovts_08_16.rds"))


    # Load occupancy data -----------------------------------------------------

    # Load site data
    sitedata <- readRDS(paste0(birdie_dir,"data/site_dat_sa_wcovts_08_19.rds"))

    # Load visit data
    visitdata <- readRDS(paste0(birdie_dir, "data/visit_dat_", sp_sel, "_wcovts_08_16.rds"))


    # Format to occuR ---------------------------------------------------------

    occuRdata <- prepDataOccuR(sitedata, visitdata)


    # Fit occupancy model -----------------------------------------------------

    # Define site model
    sitemod <- c("1", "s(water, bs = 'cs')", "s(prcp, bs = 'cs')", "s(tmax - tmin, bs = 'cs')",
                 "t2(lon, lat, occasion, k = c(15, 3), bs = c('ts', 'cs'), d = c(2, 1))")
    # The original objective was 25 knots for the spatial effect (although I ran into memory issues)

    # Define visit model
    visitmod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

    # Scale variables
    # visitvars <- visitvars %>%
    #     mutate(across(.col = -c(lon, lat, year, month, Pentad, obs, site, occasion, visit), .fns = ~scale(.x)))

    # Smooth for spatial effect on psi
    fit <- fit_occu(forms = list(reformulate(visitmod, response = "p"),
                                 reformulate(sitemod, response = "psi")),
                    visit_data = occuRdata$visit,
                    site_data = occuRdata$site)

    saveRDS(fit, paste0("/drv_birdie/birdie_FTP/", sp_sel, "/", sp_sel, "_occur_fit_08_16.rds"))

}
