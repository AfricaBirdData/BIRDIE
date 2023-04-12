library(BIRDIE)

rm(list = ls())

config <- BIRDIE::configPreambAbu(2021, server = FALSE,
                                   mod_file = "cwac_ssm_lat_season_multi_hier.R")

for(i in 1:length(config$species)){

    # Create log
    if(i == 1){
        createLog(config, log_file = NULL, date_time = NULL, species = NA, model = NA,
                  year = NA, data = NA, fit = NA, diagnose = NA, summary = NA,
                  package = NA, notes = "Log file created")
    }

    sp_code <- config$species[i]

    message(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))


    # PREPARE SSM SPECIES DATA ------------------------------------------------

    # Read in catchment data. This should go as an argument
    catchment <- sf::read_sf(file.path(config$data_dir, "quinary_catchmt_22.shp"))

    # Re-project and simplify
    catchment <- catchment %>%
        dplyr::select(QUATERNARY, Province, UNIT_ID) %>%
        sf::st_simplify(preserveTopology = TRUE, dTolerance = 1000) %>%
        sf::st_transform(crs = sf::st_crs(4326))


    # Run abudance pipeline module 1
    status_abu1 <- ppl_run_pipe_abu1(sp_code, config, steps = c("fit", "diagnose"),
                                     prep_data_steps = c("missing", "gee", "subset", "model"),
                                     catchment = catchment,
                                     force_catchm = FALSE,
                                     force_gee = FALSE,
                                     monitor = TRUE)

    message(paste("ABU1 status =", status_abu1))

}
