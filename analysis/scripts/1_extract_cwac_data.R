# 17-05-2021

# In this script we download CWAC data for the Barberspan

library(CWAC)
library(tidyverse)
devtools::load_all() # To load BIRDIE package functions

rm(list = ls())


# Extract data ------------------------------------------------------------

# We will need the migratory status of the different species later, so we need
# to load these data
bird_ids <- read.csv("analysis/data/barberspanIDs.csv")

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
counts <- map(surveys, ~.x[["records"]])

# Check that all surveys have the same fields
BIRDIE::checkListEqual(lapply(counts, names))

# If TRUE bind data frame
counts <- do.call("rbind", counts)

# Extract survey metadata
info <- map(surveys, ~.x[["summary"]])

# Check that all summaries have the same fields
BIRDIE::checkListEqual(lapply(info, names))

# Use data.table to bind dataframes with different fields
info <- data.table::rbindlist(info, fill = TRUE)

# Migratory status?? (can we find this somewhere?)
searchCwacTerm("SpeciesList")



# Combine data and info
counts <- full_join(counts, info, by = c("card" = "Card"))

# Add migratory status
counts <- left_join(counts, bird_ids, by = c("spp" = "id"))


# Save data ---------------------------------------------------------------

saveRDS(counts, "analysis/data/cwacdata.rds")

