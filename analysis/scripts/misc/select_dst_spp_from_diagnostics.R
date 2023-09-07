library(BIRDIE)
library(dplyr)

rm(list = ls())


# Configuration -----------------------------------------------------------

# Set the same configuration used for running the pipeline
# Note that models are not important, only paths to directories
det_mods <- list(det_mod1 = c("(1|obs_id)", "log_hours", "prcp", "tdiff", "cwac"),
                 det_mod2 = c("(1|obs_id)", "(1|site_id)", "log_hours", "prcp", "tdiff", "cwac"))

config <- configPipeline(year = 2010,
                         dur = 3,
                         occ_mod = c("log_dist_coast", "elev", "log_hum.km2", "wetcon",
                                     "watrec", "watext", "log_watext", "watext:watrec",
                                     "ndvi", "prcp", "tdiff"),
                         det_mod = det_mods$det_mod1,
                         fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "wetext_2018","wetcon_2018",
                                        "dist_coast", "elev"),
                         package = "spOccupancy",
                         data_dir = "analysis/hpc/imports",
                         out_dir = "analysis/hpc/imports",
                         server = TRUE)

# Select species
sp_codes <- config$species
sp_codes <- c(sp_codes, 566, 463, 220, 4127, 461, 764, 653, 626, 613, 4125, 1037, 486, 479)

# Inspect diagnostics
diags <- combineOccuDiags(config, sp_codes, 2009)

# Look for problem species
spp <- selectSppFromDiag(config, sp_codes = sp_codes, 2009)


diags %>%
    filter(sp %in% spp$no_converge)

paste0(spp$bad_fit, collapse = ", ")
