library(dplyr)

rm(list = ls())

# Set the parent directory where all model fits can be found.
# The parent directory needs to host all of the species directories
files_dir <- "/drv_birdie/birdie_ftp"

# Set the directory where the file with the model fit names is stored
out_dir <- "/home/birdie/analysis"

# Define the pattern we want to find. It can contain regular expressions
tofind <- "^occu_site_dat_"

# Extract file names without full path to handle more easily
ff <- list.files(path = files_dir,
                 pattern = tofind,
                 full.names = FALSE,
                 recursive = TRUE)

# Separate directory (species) from file name and create a data frame
if(any(grep("/", ff))){

    ff_split <- strsplit(ff, "/")

    ff_df <- data.frame(sp_code = sapply(ff_split, "[[", 1),
                        file = sapply(ff_split, "[[", 2))

    # THE SECTIONS BELOW NEED TO BE MODIFIED ACCORDING TO THE ACTUAL FILE PATHS IN
    # THE DATAFRAME AND WHAT WE WANT TO ACHIEVE.

    # In this example we want to isolate the year to be able to filter by that field
    ff_df <- ff_df %>%
        mutate(year = gsub(tofind, "", file)) %>%
        mutate(year = gsub("\\_.*\\_ZA.csv", "", year)) %>%
        mutate(year = as.integer(year)) %>%
        mutate(file_full = paste(sp_code, file, sep =  "/"))

    # Extract the last diagnostic file for each species
    keep <- ff_df %>%
        group_by(sp_code) %>%
        #filter(year == max(year)) %>%
        ungroup() %>%
        pull(file_full)


    # If we needed the full paths of the files (we tipically don't for transfers)
    ff_full <- list.files(path = files_dir,
                          pattern = tofind,
                          full.names = TRUE,
                          recursive = TRUE)


    # Select files to keep (change to ff_full if needed)
    ff_keep <- ff[ff %in% keep]

} else {

    ff_keep <- ff

}

writeLines(ff_keep, con = file.path(out_dir, "path_list.txt"))
