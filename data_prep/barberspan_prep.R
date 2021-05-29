# Load data ---------------------------------------------------------------

# Find location code
loc_code <- CWAC::listCwacSites("North West") %>%
    dplyr::filter(Name == "Barberspan") %>%
    dplyr::pull(Loc_code)

barberspan <- CWAC::getCwacSiteCounts(loc_code)


# Save data ---------------------------------------------------------------

# Save as data
usethis::use_data(barberspan)
