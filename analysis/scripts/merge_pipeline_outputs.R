library(BIRDIE)

rm(list = ls())

modules <- c("abu", "dst")


# Merge abundance estimates -----------------------------------------------

if("abu" %in% modules){

    # Configure
    config <- configPipeline(
        year = 2021,
        dur = 29,
        mod_file = "cwac_ssm_two_season_mean_rev.R",
        package = "jagsUI",
        data_dir = NULL,
        out_dir = NULL,
        server = TRUE
    )

    # Merge and save in exports
    createCombinedExportFile(config, type = "abu")

}



# Merge occupancy estimates -----------------------------------------------

if("dst" %in% modules){

    # Configure
    config <- configPipeline(
        year = 2019,
        dur = 12,
        package = "spOccupancy",
        data_dir = NULL,
        out_dir = NULL,
        server = TRUE
    )

    # Merge and save in exports
    createCombinedExportFile(config, type = "dst")

}
