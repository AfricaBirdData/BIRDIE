library(rredlist)
library(tidyverse)

help("rredlist")

countries <- rl_countries(key = Sys.getenv("IUCN_REDLIST_KEY"))[[2]]

countries[countries$country == "South Africa",]

sa_spp <- rl_sp_country(country = "ZA", key = Sys.getenv("IUCN_REDLIST_KEY"))$result

bbp_spp <- BIRDIE::barberspan

bbp_iucn <- bbp_spp %>%
    mutate(scientific_name = paste(taxon.Genus, taxon.Species)) %>%
    dplyr::select(scientific_name, spp) %>%
    left_join(dplyr::select(sa_spp, scientific_name, category))

# Unfortunately scientific names are not exactly the same so there are errors

summary(bbp_iucn)

missed <- is.na(bbp_iucn$category)

bbp_iucn[missed,]

# Try fuzzy matching
bbp_iucn <- fuzzyjoin::stringdist_join(
    bbp_spp %>%
        mutate(scientific_name = paste(taxon.Genus, taxon.Species)) %>%
        dplyr::select(scientific_name, spp),
    sa_spp %>%
        dplyr::select(scientific_name, category),
    by = "scientific_name",
    mode = "left",
    ignore_case = FALSE,
    method = "jw",
    max_dist = 1,
    distance_col = "dist") %>%
    group_by(scientific_name.x) %>%
    slice_min(order_by = dist, n = 1)

bbp_iucn %>%
    distinct(scientific_name.x, scientific_name.y, .keep_all = TRUE) %>%
    # filter(dist > 0.1) %>%
    print(n = Inf)

bbp_iucn %>%

sa_spp[grep("cristata$",sa_spp$scientific_name), ]

# WE PROBABLY SHOULDN'T USE FUZZY MATCHING WITH GREATER DISTANCE THAN 0.09


bbp_iucn <- bbp_iucn %>%
    distinct(scientific_name.x, spp, scientific_name.y, .keep_all = TRUE) %>%
    mutate(category = if_else(dist > 0.09, NA_character_, category))

bbp_iucn %>%
    print(n = Inf)

bbp_iucn %>%
    group_by(category) %>%
    count(category)

bbp_iucn %>%
    filter(category == "NT")
