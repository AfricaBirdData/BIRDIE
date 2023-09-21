library(BIRDIE)

rm(list = ls())


# Configuration -----------------------------------------------------------

config <- configPipeline(
    year = 2022,
    dur = 30,
    region = "southafrica",
    module = "abu",
    mod_file = "cwac_ssm_two_season_mean_rev_jump.R",
    package = "jagsUI",
    data_dir = NULL,     # this might have to be adapted?
    out_dir = NULL,     # this might have to be adapted?
    server = FALSE
)


# Read in catchment data. This should go as an argument
catchment <- sf::read_sf(file.path(config$data_dir, "quinary_catchmt_22.shp"))

# Re-project and simplify
catchment <- catchment %>%
    dplyr::select(QUATERNARY, Province, UNIT_ID) %>%
    sf::st_simplify(preserveTopology = TRUE, dTolerance = 1000) %>%
    sf::st_transform(crs = sf::st_crs(4326))


# Create logs -------------------------------------------------------------

# Create log?
createLog(config, log_file = NULL, date_time = NULL, species = NA, model = NA,
          year = NA, data = NA, fit = NA, diagnose = NA, summary = NA,
          package = NA, notes = "Log file created")


# Run modules -------------------------------------------------------------

for(i in 1:length(config$species)){

    sp_code <- config$species[i]

    message(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

    # Run abudance pipeline module 1
    status_abu1 <- ppl_run_pipe_abu1(sp_code, config,
                                     steps = c("data", "fit", "diagnose", "summary"),
                                     prep_data_steps = c("subset", "missing", "gee", "model"),
                                     summary_scale = "linear",
                                     catchment = catchment,
                                     force_gee_upload = FALSE,
                                     force_gee = FALSE,
                                     monitor = TRUE)

    message(paste("ABU1 status =", status_abu1))

}
