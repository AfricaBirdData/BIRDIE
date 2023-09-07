# With this script we create a file to store indicator values in all species directories

library(dplyr)

rm(list = ls())

# Should we remove existing files first
overwrite <- TRUE

# Data directory
data_dir <- "analysis/out_nosync"

# Species directories
spp_dirs <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)[1:5]

# indicator files
indtr_files <- file.path(data_dir, spp_dirs, paste0("indtr_", spp_dirs, ".csv"))

if(overwrite){
    file.remove(indtr_files)
}

indtr_df <- data.frame(species = character(),
                       indicator = character(),
                       start_date = character(),
                       end_date = character(),
                       estimate = numeric(),
                       st_dev = numeric(),
                       lb95 = numeric(),
                       ub95 = numeric())

# Species indicator files
indtr_files <- indtr_files[!file.exists(indtr_files)]

# Create files if they don't exist
for(i in seq_along(indtr_files)){
    write.csv(indtr_df, indtr_files[i], row.names = FALSE)
}
