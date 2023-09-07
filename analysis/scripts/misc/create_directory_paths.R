rm(list = ls())

# Species
sp_codes <- unique(BIRDIE::waterbirds$SppRef)
sp_codes <- sp_codes[!is.na(sp_codes)]
sp_codes <- c(sp_codes, 566, 463, 220, 4127, 461, 764, 653, 626, 613, 4125, 1037, 486, 479)

# Define source and target directories
out_dir <- ""
in_dir <- "analysis/hpc/imports"
# in_dir <- "/drv_birdie/birdie_ftp"

# File name
prefix <- "occu_diagnostics_"
suffix <- "_ZA.csv"
year_sel <- "2019"
# package <- "jagsUI"
package <- "spOccupancy"

# Create model output import-from file names -------------------------------

# ff_out <- vector(length = length(sp_codes)*length(config$years))
ff_out <- vector(length = length(sp_codes))
i <- 1
for(sp_code in sp_codes){
    filename <- paste0(prefix, package, "_", year_sel, "_", sp_code, suffix)
    ff_out[i] <- file.path(out_dir, sp_code, filename)
    # ff_out[i] <- file.path(out_dir, sp_code)
    i <- i+1
    # }
}

writeLines(ff_out, con = "analysis/out_nosync/output_filenames.txt")
