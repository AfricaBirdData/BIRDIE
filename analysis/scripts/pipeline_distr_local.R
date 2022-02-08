library(BIRDIE)
library(dplyr)
# library(occuR)

rm(list = ls())

test_years <- c(2009, 2017, 2018)
config <- configPreambOccuR(year = 2008, server = FALSE)

# Create indicator storage files?
ppl_create_indtr_files(config, overwrite_indtr = TRUE)


# Overwrite?
overwrite_indtr <- FALSE

# Create indicator file path
indtr_path <- file.path(config$fit_dir, "group", paste0("indtr_", "group", ".csv"))

# Check if file exists
if(!file.exists(indtr_path) | (file.exists(indtr_path) & overwrite_indtr)){

    # Create empty data frame
    indtr <- data.frame(species = character(),
                        indicator = character(),
                        start_date = character(),
                        end_date = character(),
                        term = character(),
                        estimate = numeric(),
                        st_dev = numeric(),
                        lb95 = numeric(),
                        ub95 = numeric(),
                        opt = numeric())

    # Save
    write.csv(indtr, indtr_path, row.names = FALSE)

}

for(y in seq_along(test_years)){

    year <- test_years[y]
    config <- configPreambOccuR(year = year, server = FALSE)

    for(i in seq_along(config$species)){

        sp_code <- config$species[i]

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
                           print_fitting = TRUE,
                           verbose = TRUE)

    }

    # Estimate species type indicators
    ppl_estimate_distr_sp_type(sp_type = "group", year, config, force_predict = FALSE, verbose = TRUE)

}
