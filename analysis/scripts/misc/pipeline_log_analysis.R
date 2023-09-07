library(dplyr)
library(lubridate)
library(ggplot2)
library(BIRDIE)


rm(list = ls())

# Read data in
logs <- read.csv("analysis/hpc/reports/pipe_log_dst_2023-04-24.csv")


# Convert dates to date-time objects and filter out anything that are not model fits
logs <- logs %>%
    mutate(date_time = as_datetime(date_time)) %>%
    filter(!is.na(fit))

logs %>%
    count(species) %>%
    filter(n < 2)

logs %>%
    filter(fit != 0) %>%
    arrange(species, year)


# Find duration of model fit
logs <- logs %>%
    arrange(species, date_time) %>%
    group_by(species) %>%
    mutate(dur = difftime(date_time, lag(date_time), units = "hours")) %>%
    ungroup()


logs %>%
    filter(dur < 30) %>%
    ggplot() +
    geom_histogram(aes(x = dur))

logs %>%
    group_by(species) %>%
    summarise(total_dur = sum(dur, na.rm = TRUE)) %>%
    ggplot() +
    geom_histogram(aes(x = total_dur))


logs %>%
    filter(dur < 30) %>%
    pull(dur) %>%
    as.numeric() %>%
    summary()

logs %>%
    filter(dur < 30, dur > 8)

# Define species to fit models to
total_spp <- unique(BIRDIE::barberspan$SppRef) # For now, we want to select species present at Barberspan

# Remove partially identified species
total_spp <- total_spp[total_spp < 10000]

logs %>%
    filter(year == 2008) %>%
    arrange(species) %>%
    pull(species) %>% unique() %>% length()

length(unique(barberspan$Species))

# Analyse time excess reports ---------------------------------------------

ff <- list.files("analysis/hpc/reports/reports_2023_04_24")

# subset computation time ones
ff <- ff[grep("Computation_time", ff)]

# Extract species and years
spp <- gsub("Computation_time_exceeded_sp_", "", ff)
spp <- gsub("_2.*.txt", "", spp)

yy <- gsub("Computation_time_exceeded_sp_.*_", "", ff)
yy <- gsub(".txt", "", yy)

# Pair species and years
to_run <- data.frame(sp = as.integer(spp), year = as.integer(yy))

to_run %>%
    filter(year == 2010) %>%
    arrange(sp) %>%
    pull(sp) %>% length() paste(collapse = ",")
