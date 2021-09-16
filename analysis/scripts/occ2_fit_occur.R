library(occuR)
library(tidyverse)
library(BIRDIE)

rm(list = ls())


# For now, we want to select species present at Barberspan
bbpan <- BIRDIE::barberspan %>%
    pull(spp) %>%
    unique()


# Select a species and a region -------------------------------------------

# We are most interested in the Maccoa Duck VU (103) and the Cape Cormorant EN (48)
# i <- 1
# sp_sel <- bbpan[i]
sp_sel <- 103

sp_name <- BIRDIE::barberspan %>%
    dplyr::filter(spp == sp_sel) %>%
    mutate(name = paste(taxon.Common_species, taxon.Common_group)) %>%
    pull(name) %>%
    unique()


# Load occupancy data -----------------------------------------------------

# Load site data
sitedata <- readRDS("analysis/data/site_dat_sa_wcovts_16_19.rds")

# Load visit data
visitdata <- readRDS(paste0("analysis/data/visit_dat_", sp_sel, "_wcovts_16_19.rds"))


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

saveRDS(fit, paste0("analysis/output/", sp_sel, "_occur_fit_16_19.rds"))


