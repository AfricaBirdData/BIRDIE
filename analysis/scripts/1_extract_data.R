# 17-05-2021

# In this script we download CWAC data for the Barberspan

library(CWAC)
library(tidyverse)

rm(list = ls())


# Extract data ------------------------------------------------------------

# Find location code
loc_code <- listCwacSites("North West") %>%
  filter(Name == "Barberspan") %>%
  pull(Loc_code)

# Extract location cards
cards <- listCwacCards(loc_code)

# Download surveys based on location cards
surveys <- map(cards$Card, ~getCwacSurvey(.x))


# Join counts with survey metadata ----------------------------------------

# Extract counts
counts <- map_df(surveys, ~.x[["records"]])

# Extract survey metadata
info <- map_df(surveys, ~.x[["summary"]])



# Migratory status?? (can we find this somewhere?)
searchCwacTerm("SpeciesList")



# Combine data and info
counts <- left_join(counts, info, by = c("card" = "Card"))



# Save data ---------------------------------------------------------------

saveRDS(counts, "analysis/data/barberspan.rds")

