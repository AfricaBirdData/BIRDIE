library(dplyr)
library(readxl)

rm(list = ls())

# CWAC list to compare with
# spp <- CWAC::listCwacSpp()

waterbirds <- read_excel("analysis/data/BIRDIE Project Waterbird master list_15052023.xlsx",
                      sheet = "Masterlist")

# Correct variable names and make them match CWAC database
waterbirds <- waterbirds %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    rename(SppRef = `#_Ref`) %>%
    select(SppRef, everything())


# Save data ---------------------------------------------------------------

# Save as data
usethis::use_data(waterbirds, overwrite = TRUE)
