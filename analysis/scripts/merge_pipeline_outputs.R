library(BIRDIE)

rm(list = ls())

modules <- c("abu", "dst", "env")


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


# Merge environmental data ------------------------------------------------

if("env" %in% modules){

    # Configure
    config <- configPipeline(
        year = 2022,
        dur = 14,
        package = "spOccupancy",
        data_dir = NULL,
        out_dir = NULL,
        server = TRUE
    )

    # Model period (years)
    years_per_file <- 3

    periods <- with(config, seq(year - dur, year, years_per_file))
    periods <- paste(periods, periods+years_per_file-1, sep = "_")
    periods <- gsub("^.{0,2}|\\_.{0,2}", "_", periods)

    # Merge and save in exports
    fpaths <- file.path(config$out_dir, paste0("site_dat_sa_gee", periods, ".csv"))

    env_out <- data.frame()

    for(i in 1:length(fpaths)){

        env_file <- fpaths[i]

        if(file.exists(env_file)){

            new_env <- utils::read.csv(env_file) %>%
                dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3)))

            env_out <- dplyr::bind_rows(env_out, new_env)

        } else {
            env_out <- env_out
        }

    }

    utils::write.csv(env_out,
                     setSpOutFilePath("site_dat_gee", config, config$years_ch, "export", ".csv"),
                     row.names = FALSE)

}
