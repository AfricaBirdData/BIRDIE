library(dplyr)

rm(list = ls())

# Set the parent directory where all model fits can be found.
# The parent directory needs to host all of the species directories
files_dir <- "analysis/out_nosync"

# Extract file names without full path to handle more easily
ff <- list.files(path = files_dir,
                 pattern = "^occu_fit",
                 full.names = FALSE,
                 recursive = TRUE)

# Separate directory (species) from file name and create a data frame
ff_split <- strsplit(ff, "/")


ff_df <- data.frame(sp_code = sapply(ff_split, "[[", 1),
                    file = sapply(ff_split, "[[", 2)) %>%
    mutate(year = gsub("occu_fit_spOccupancy_", "", file)) %>%
    mutate(year = gsub("\\_.*", "", year)) %>%
    mutate(year = as.integer(year)) %>%
    mutate(file_full = paste(sp_code, file, sep =  "/"))

# Extract the last model fit per species
keep <- ff_df %>%
    group_by(sp_code) %>%
    filter(year == max(year)) %>%
    ungroup() %>%
    pull(file_full)


# Now we need the full paths of the files
ff_full <- list.files(path = files_dir,
                      pattern = "^occu_fit",
                      full.names = TRUE,
                      recursive = TRUE)


# Select files to keep
ff_keep <- ff_full[ff %in% keep]

writeLines(ff_keep, con = "analysis/out_nosync/files_to_keep.txt")
