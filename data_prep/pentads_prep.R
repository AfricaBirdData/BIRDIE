library(tidyverse)
library(sf)

rm(list = ls())

# Load pentad data from kml. This was downloaded from the website
pentads <- sf::st_read("../resources/project_sabap2.kml")

# Remove description because it is empty
pentads <- dplyr::select(pentads, -Description)

# Remove Z dimension
pentads <- sf::st_zm(pentads)

# Save as data
usethis::use_data(pentads)


# Modify existing data ----------------------------------------------------

pentads <- BIRDIE::pentads

# Remove Z dimension
pentads <- sf::st_zm(pentads)

usethis::use_data(pentads, overwrite = TRUE)
