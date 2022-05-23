
# In this script we take Doug's Du Toits Pan waterbird data and format it as
# CWAC data to feed into the pipeline

library(dplyr)

rm(list = ls())


# Download CWAC data ------------------------------------------------------

# # Download Du Toits Pan CWAC data and save if necessary
# counts <- CWAC::getCwacSiteCounts(28462448)
# saveRDS(counts, "analysis/dutoitspan_data/dutoitspan_cwac_data.rds")

# # Get CWAC species list and save
# cwac_spp <- CWAC::listCwacSpp()
# saveRDS(cwac_spp, "analysis/dutoitspan_data/cwac_spp_list.rds")


# Read in data ------------------------------------------------------------

# Read in Dougs data
dougs <- readr::read_csv("analysis/data/dutoitspan_data/Waterbird data_Doug.csv")

# Read in CWAC data
cwac_data <- readRDS("analysis/data/dutoitspan_data/dutoitspan_cwac_data.rds")

# Read in CWAC species list
cwac_spp <- readRDS("analysis/data/dutoitspan_data/cwac_spp_list.rds")


# Clean data --------------------------------------------------------------

# Remove extra variables (2019/10/01 was the last survey)
which(names(dougs) == "2019/10/01")

dougs <- dougs[, 1:34]

summary(dougs)

# The last row seems to be NA

# How many NA values in each row?
apply(dougs, 1, function(x) sum(is.na(x)))

# Remove NA rows
dougs <- dougs[complete.cases(dougs),]


# Add species code to Doug's counts ---------------------------------------

# Add species code to Doug's data
dougs <- dougs %>%
    left_join(cwac_spp %>%
                  mutate(Common_group = if_else(is.na(Common_group), "", Common_group),
                         Species = paste(Common_species, Common_group),
                         Species = gsub(" $", "", Species)) %>%
                  select(Species, SppRef),
              by = "Species")

dougs[is.na(dougs$SppRef),]

cwac_spp[grep("Teal", cwac_spp$Common_group),]
cwac_spp[grep("Ibis", cwac_spp$Common_group),]

# Change species names to match those of CWAC
dougs <- dougs %>%
    mutate(Species = case_when(Species == "Hottentot Teal" ~ "Blue-billed Teal",
                               Species == "Hadeda Ibis" ~ "Hadada Ibis",
                               TRUE ~ Species))

# Re-join species with species codes
dougs <- dougs %>%
    select(-SppRef) %>%
    left_join(cwac_spp %>%
                  mutate(Common_group = if_else(is.na(Common_group), "", Common_group),
                         Species = paste(Common_species, Common_group),
                         Species = gsub(" $", "", Species)) %>%
                  select(Species, SppRef),
              by = "Species") %>%
    mutate(SppRef = as.integer(SppRef))

dougs[is.na(dougs$SppRef),]

# Add Common_group and Common_species columns to Doug's data
nn <- strsplit(dougs$Species, " ")

dougs <- dougs %>%
    mutate(Common_group = sapply(nn, function(x) x[length(x)]),
           Common_species = sapply(nn, function(x) paste(x[1:(length(x)-1)], collapse = " ")))


# Transform to long format and fix dates ----------------------------------

# Convert to long format
dougs_l <- dougs %>%
    dplyr::select(Species, SppRef, everything()) %>%
    tidyr::pivot_longer(cols = -c(Species, SppRef, Common_group, Common_species), names_to = "StartDate", values_to = "Count")

# Transform dates
dougs_l <- dougs_l %>%
    mutate(#StartDate = gsub("X", "", StartDate),
           StartDate = as.Date(StartDate, format = "%Y/%m/%d"))

# Add seasons
dougs_l <- dougs_l %>%
    mutate(Year = lubridate::year(StartDate),
           month = lubridate::month(StartDate),
           day = lubridate::day(StartDate)) %>%
    mutate(Season = case_when((month == 1) & day > 14 ~ "S",
                              (month == 2) & day < 16 ~ "S",
                              month == 7 ~ "W",
                              TRUE ~ "O")) %>%
    select(-c(month, day))


# Bind CWAC and Doug's data -----------------------------------------------

final_counts <- bind_rows(cwac_data[0,], dougs_l)

summary(final_counts)

# Fix X and Y columns because they will be needed later
final_counts$X <- cwac_data$X[1]
final_counts$Y <- cwac_data$Y[1]

# Add LocationCode and name
final_counts$LocationCode <- unique(cwac_data$LocationCode)
final_counts$LocationName <- unique(cwac_data$LocationName)
final_counts$Province <- unique(cwac_data$Province)
final_counts$Country <- unique(cwac_data$Country)

# Save
write.csv(final_counts, "analysis/data/28462448_data_2022_doug.csv", row.names = FALSE)
