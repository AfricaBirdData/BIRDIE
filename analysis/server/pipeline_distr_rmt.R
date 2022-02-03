library(BIRDIE)
library(dplyr)
# library(occuR)

rm(list = ls())

test_years <- c(2010, 2017)
config <- configPreambOccuR(year = 2017, server = TRUE)

# Create indicator storage files?
ppl_create_indtr_files(config, overwrite_indtr = TRUE)

for(i in seq_along(config$species)){
    for(y in seq_along(test_years)){

        sp_code <- config$species[i]
        year <- test_years[y]

        config <- configPreambOccuR(year = year, server = TRUE)

        # Species name
        sp_name <- BIRDIE::barberspan %>%
            dplyr::filter(SppRef == sp_code) %>%
            mutate(name = paste(Common_species, Common_group)) %>%
            pull(name) %>%
            unique()

        print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

        ppl_run_pipe_distr(sp_code = sp_code,
                           sp_name = sp_name,
                           year = year,
                           config = config,
                           steps = c("indtr"),
                           download_from_abap = TRUE,
                           save_occu_data = TRUE,
                           overwrite_occu_data = c("site", "visit", "det"),
                           scale_vars_occur = list(visit = NULL,
                                                   site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")),
                           print_fitting = FALSE,
                           verbose = TRUE)

    }
}
