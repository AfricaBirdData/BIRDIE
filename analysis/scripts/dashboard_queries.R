library(BIRDIE)
library(dplyr)
library(ggplot2)

rm(list = ls())

test_years <- c(2009, 2017, 2018)
config <- configPreambOccuR(year = 2008, server = FALSE)

# remove 240 for now
config$species <- c(4, 6, 41, 235)

# indicators <- list(distr = c("occu", "aoo", "daoo", "abu", "dabu"))
# species <- c(4, 6, 41, 235, 240)
# term <- c("annual", "short", "long")

indtr_orig <- read.csv("analysis/scripts/indicators/diagrams/indtr_all.csv")


# Refine queries ----------------------------------------------------------

# If at site only show over time because users can hover over plot
indtr <- indtr_orig %>%
    filter(cond1 == "across" | cond2 == "across" | cond3 == "across")

indtr


# Create group indicator directory ----------------------------------------

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


# Estimate species type occupancy -----------------------------------------


for(y in seq_along(test_years)){

    year <- test_years[y]

    config <- configPreambOccuR(year = year, server = FALSE)

    ppl_estimate_distr_sp_type(sp_type = "group", year, config, force_predict = FALSE, verbose = TRUE)

}


# Test function -----------------------------------------------------------

source("beta/generateQuery.R")

generateQuery(.target = "species", .sites = "national", .years = 2010, .indtr = "occu", .term = "annual", .species = 4)

generateQuery(.target = "species", .species = 4, .sites = "national", .years = c(2008, 2010), .indtr = "aoo", .term = "annual")

generateQuery(.target = "species", .species = 4, .sites = "national", .years = 2010, .indtr = "daoo", .term = "annual")

generateQuery(.target = "species", .species = 4, .sites = 26352535, .years = 2010, .indtr = "abu", .term = "annual")

generateQuery(.target = "species_type", .sites = "national", .years = 2010, .indtr = "occu", .term = "annual", .sp_type = "group", config = config)

generateQuery(.target = "species_type", .sites = "national", .years = 2009, .indtr = "occu", .term = "annual", .sp_type = "group", config = config)

?ppl_estimate_aoo
ppl_estimate_aoo("group", 2009, config, verbose = TRUE)

# Remove last species from config
config$species <- species <- c(4, 6, 41, 235)

source("beta/estimateSpeciesTypeIndtr.R")
estimateSpeciesTypeIndtr(sp_type = "group", year = 2009, config = config, verbose = TRUE)
estimateSpeciesTypeIndtr(sp_type = "group", year = 2019, config = config, verbose = TRUE)


ppl_estimate_daoo("group", 2019, config, term = "short", verbose = TRUE)
ppl_estimate_daoo("group", 2019, config, term = "annual", verbose = TRUE)
