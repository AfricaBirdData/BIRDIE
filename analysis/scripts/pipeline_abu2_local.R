library(BIRDIE)
library(dplyr)

rm(list = ls())

year <- 2017
config <- configPreambJAGS(year = year, server = FALSE)
sites <- CWAC::listCwacSites(.region_type = "country", .region = "South Africa")

site_ids <- sites %>%
    pull(LocationCode) %>%
    unique()

# for(s in seq_along(site_ids)){

    # site_id <- site_ids[s]
    site_id <- 26352535

    # Run branch 1 of abundance pipeline for site
    ppl_run_pipe_abu2(site = site_id,
                      year = year,
                      config = config,
                      steps = c("fit", "summ"))


# }
