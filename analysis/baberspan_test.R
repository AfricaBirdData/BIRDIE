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

# Download cards
counts <- map_df(cards$Card, ~getCwacSurvey(.x)[["records"]])

# Download metadata
info <- map_df(cards$Card, ~getCwacSurvey(.x)[["summary"]])

# Combine data and info
counts <- left_join(counts, info, by = c("card" = "Card"))

# Save data
saveRDS(counts, "analysis/data/barberspan.rds")

