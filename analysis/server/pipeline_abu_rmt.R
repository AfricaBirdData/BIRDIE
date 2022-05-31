library(BIRDIE)

rm(list = ls())

config <- BIRDIE::configPreambJAGS(2017, server = TRUE)

for(i in 1:length(config$species)){

    sp_code <- config$species[i]

    print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))


    # PREPARE SSM SPECIES DATA ------------------------------------------------

    # Read in catchment data. This should go as an argument
    catchment <- sf::read_sf(file.path(config$data_dir, "catchmt_4.shp"))

    # Remove marine area around South Africa and also neighboring countries
    # and simplify
    catchment <- catchment %>%
        dplyr::filter(AREA < 150) %>%
        dplyr::select(QUATERNARY, QUAT_CODE) %>%
        sf::st_simplify(preserveTopology = TRUE, dTolerance = 1000)

    counts <- ppl_create_data_ssm(sp_code, config$year, catchment, config,
                                  steps = c("missing", "gee", "subset"),
                                  upload_catchment = FALSE, force_gee = TRUE)


    # Fit model ---------------------------------------------------------------

    ppl_fit_ssm_model(sp_code, config)

    ppl_summarize_ssm(sp_code, config)

}
