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

# Select a range of time. Occupancy models will be fitted from first to second
ini_year <- 2008
years <- c(ini_year, ini_year + 5)
years_ch <- paste(substring(as.character(years), 3, 4), collapse = "_")

# for(i in seq_along(bbpan)){
i = 1

    # Select a species and a region -------------------------------------------

    sp_sel <- bbpan[i]

    print(paste0("Working on species ", sp_sel, " (", i, " of ", length(bbpan), ")"))

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
                                  years = years[1]:years[2],
                                  clim_covts = c("prcp", "tmax", "tmin", "aet", "pet"),
                                  covts_dir = paste0(birdie_dir, "downloads/"),
                                  file_fix = c("terraClim_", "_03_19"),
                                  savedir = paste0(birdie_dir, "data/pentads_sa.rds"))

    future::plan("sequential")

    saveRDS(visitdata, paste0(birdie_dir, "data/visit_dat_", sp_sel, "_wcovts_", years_ch, ".rds"))


    # Format to occuR ---------------------------------------------------------

    occuRdata <- prepDataOccuR(sitedata, visitdata)


    # Fit occupancy model -----------------------------------------------------

    # Define site model
    sitemod <- c("1", "s(water, bs = 'cs')", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')",
                 "t2(lon, lat, occasion, bs = c('ts', 'cs'), d = c(2, 1))")

    # Define visit model
    visitmod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

    # Scale variables
    occuRdata$visit <- occuRdata$visit %>%
        mutate(tdiff = tmax - tmin,
               across(.col = c(prcp, tdiff), .fns = ~scale(.x)))

    print(paste0("Fitting model at ", Sys.time(), ". This will take a while..."))

    # Smooth for spatial effect on psi
    fit <- fit_occu(forms = list(reformulate(visitmod, response = "p"),
                                 reformulate(sitemod, response = "psi")),
                    visit_data = occuRdata$visit,
                    site_data = occuRdata$site,
                    print = TRUE)

    saveRDS(fit, paste0("/drv_birdie/birdie_ftp/", sp_sel, "/", sp_sel, "_occur_fit_", years_ch, ".rds"))

# }
