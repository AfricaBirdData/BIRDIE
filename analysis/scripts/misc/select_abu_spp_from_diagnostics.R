library(BIRDIE)
library(dplyr)

rm(list = ls())


# Configuration -----------------------------------------------------------

# Set the same configuration used for running the pipeline
# Note that models are not important, only paths to directories
config <- configPipeline(
    year = 2021,
    dur = 29,
    mod_file = "cwac_ssm_two_season_mean_rev.R",
    package = "jagsUI",
    data_dir = NULL,
    out_dir = "analysis/downloads",
    server = TRUE
)

# Select problem species
sp_codes <- config$species
sp_codes <- c(566, 463, 220, 4127, 461, 764, 653, 626, 613, 4125, 1037, 486, 479)

# Inspect diagnostics
diags <- combineAbuDiags(config, sp_codes, config$years_ch)

diags <- diags %>%
    mutate(nc_rate = nc_obs/nobs)

diags %>%
    arrange(desc(nc_rate))

spp <- selectSppFromDiag(config, sp_codes = sp_codes, config$years_ch, module = "abu")
