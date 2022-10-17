library(BIRDIE)

rm(list = ls())

config <- BIRDIE::configPreambJAGS(2017, server = TRUE)

for(i in 1:length(config$species)){

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


    status_abu1 <- ppl_run_pipe_abu1(sp_code, config, steps = c("data", "fit", "summary"),
                                     prep_data_steps = c("missing", "gee", "subset", "model"),
                                     catchment = catchment,
                                     upload_catchment = FALSE, force_gee = TRUE)

    message(paste("ABU1 status =", status_abu1))

}
