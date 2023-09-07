library(BIRDIE)

# Load data ---------------------------------------------------------------

# Find location code
loc_code <- CWAC::listCwacSites("North West") %>%
    dplyr::filter(LocationName == "Barberspan") %>%
    dplyr::pull(LocationCode)

barberspan <- CWAC::getCwacSiteCounts(loc_code)


# Save data ---------------------------------------------------------------

# Save as data
usethis::use_data(barberspan, overwrite = TRUE)
