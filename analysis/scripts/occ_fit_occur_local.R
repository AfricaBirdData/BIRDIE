# In this script we fit an occupancy model with occuR for a group of species

library(occuR)
library(dplyr)
library(tidyr)
library(BIRDIE)
library(sf)

rm(list = ls())


# Script parameters -------------------------------------------------------

year_sel <- 2008

config <- configPreambOccuR(year = year_sel, server = FALSE)


# Load site and visit data ------------------------------------------------

# Load data and subset years
sitedata <- readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    dplyr::select(Name, lon, lat, watocc_ever, dist_coast, ends_with(match = as.character(config$years))) %>%
    tidyr::drop_na()   # I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE

visitdata <- readRDS(file.path(config$data_dir, "visit_dat_sa_gee_08_19.rds")) %>%
    filter(year %in% config$years)


# Define models -----------------------------------------------------------

# Detection
visit_mod <- c("1", "log(TotalHours+1)", "s(month, bs = 'cs')")

# Occupancy
site_mods <- list(mod1 = c("-1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec", config$sp_temp),
                  mod2 = c("1", "dist_coast", "s(prcp, bs = 'cs')", "s(tdiff, bs = 'cs')", "s(ndvi, bs = 'cs')", "s(watext, bs = 'cs')", "s(watrec, bs = 'cs')"),
                  mod3 = c("1", "dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec"))


# Fit models --------------------------------------------------------------

for(i in seq_along(config$species)){

    # This will be used for memory cleaning
    keep <- ls()

    # Select a species and a region -------------------------------------------

    sp_sel <- config$species[i]

    print(paste0("Working on species ", sp_sel, " (", i, " of ", length(config$species), ")"))


    # Format to occuR ---------------------------------------------------------

    occuRdata <- BIRDIE::prepDataOccuR(spp_code = sp_sel,
                                       years = config$years,
                                       site_data = sitedata,
                                       visit_data = visitdata,
                                       download = TRUE)


    # Fit occupancy model -----------------------------------------------------

    print(paste0("Fitting model at ", Sys.time(), ". This will take a while..."))

    # Fit models sequentially if they don't work
    success <- FALSE
    m <- 0
    while(!success && m <= length(site_mods)){
        m <- m + 1
        site_mod <- site_mods[[m]]
        print(paste("Trying model", m))
        tryCatch({
            fit <- fit_occu(forms = list(reformulate(visit_mod, response = "p"),
                                         reformulate(site_mod, response = "psi")),
                            visit_data = occuRdata$visit,
                            site_data = occuRdata$site,
                            print = TRUE)

            success <- TRUE
            saveRDS(fit, file.path(config$fit_dir, sp_sel, paste0("occur_fit_", config$years_ch, "_", sp_sel, ".rds")))
        }, error = function(e){
            success <- FALSE
            sink(file.path(config$fit_dir, sp_sel, paste0("failed_fit_", m, "_", sp_sel,".txt")))
            print(e)
            sink()}) # TryCatch fit
    }

    # Clean
    rm(list = setdiff(ls(), keep))
    gc()

}
