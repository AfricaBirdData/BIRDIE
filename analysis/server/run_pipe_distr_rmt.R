library(BIRDIE)
library(occuR)
library(dplyr)

rm(list = ls())


year = 2017
config <- configPreambOccuR(year = year, server = TRUE)
beta_dir <- "/drv_birdie/Working/git/BIRDIE"

for(i in seq_along(config$species)){

    sp_code = config$species[i]

    # Species name
    sp_name <- BIRDIE::barberspan %>%
        dplyr::filter(SppRef == sp_code) %>%
        mutate(name = paste(Common_species, Common_group)) %>%
        pull(name) %>%
        unique()


    source(file.path(beta_dir, "beta/create_site_visit.R"))
    source(file.path(beta_dir, "beta/fit_occur_model.R"))
    source(file.path(beta_dir, "beta/summarize_occur.R"))
    # source("beta/predict_occur.R")
    # source("beta/estimate_aoo.R")

    run_pipe_distr <- function(sp_code, year, ...){

        create_site_visit(sp_code, year, ...)

        fit_occur_model(sp_code, year, ...)

        summarize_occur(sp_code, year, ...)

        # estimate_aoo(sp_code, year, ...)

    }

    run_pipe_distr(sp_code = sp_code,
                   year = year,
                   config = config,
                   download_from_abap = TRUE,
                   save_occu_data = TRUE,
                   overwrite_occu_data = c("site", "visit", "det"),
                   scale_vars_occur = list(visit = NULL,
                                           site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")),
                   print_fitting = FALSE,
                   sp_name = sp_name,
                   verbose = TRUE)
}

