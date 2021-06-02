library(tidyverse)
library(sf)

rm(list = ls())

# Load pentad data from kml. This was downloaded from the website
pentads <- sf::st_read("../resources/project_sabap2.kml")

# Remove description because it is empty
pentads <- dplyr::select(pentads, -Description)

# Save as data
usethis::use_data(pentads)
