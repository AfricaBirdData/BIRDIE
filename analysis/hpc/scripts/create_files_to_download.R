library(BIRDIE)

rm(list = ls())

# Create config object
config_year <- 2010

config <- configPipeline(year = config_year,
                         dur = 29,
                         mod_file = "cwac_ssm_two_season_mean_rev.R",
                         server = TRUE,
                         package = "jagsUI",
                         data_dir = "/home/crvfra001/birdie/data",
                         out_dir = "/scratch/crvfra001/birdie/output")

sp_codes <- c(config$species, 566, 463, 220, 4127, 461, 764, 653, 626, 613, 4125, 1037, 486, 479)

writeLines(as.character(sp_codes), con = "analysis/hpc/scripts/spp_codes.txt")

# Define source and target directories
out_dir <- "/scratch/crvfra001/birdie/output"
in_dir <- "analysis/hpc/imports"
in_dir <- "/drv_birdie/birdie_ftp"
prefix <- NULL
suffix <- NULL

# Create model output import-from file names -------------------------------

# ff_out <- vector(length = length(sp_codes)*length(config$years))
ff_out <- vector(length = length(sp_codes))
i <- 1
for(sp_code in sp_codes){
    # for(t in seq_along(config$years)){
        year_sel <- config$years[t]
        filename <- paste0("occu_fit_", config$package, "_", year_sel, "_", sp_code, ".rds")
        # ff_out[i] <- file.path(out_dir, sp_code, filename)
        ff_out[i] <- file.path(out_dir, sp_code)
        i <- i+1
    # }
}

writeLines(ff_out, con = "analysis/hpc/scripts/dst_model_output_filenames.txt")



# Create model output export-to file names --------------------------------

ff_in <- vector(length = length(sp_codes)*length(config$years))
i <- 1
for(sp_code in sp_codes){
    for(t in seq_along(config$years)){
        year_sel <- config$years[t]
        filename <- paste0("occu_fit_", config$package, "_", year_sel, "_", sp_code, ".rds")
        ff_in[i] <- file.path(in_dir, sp_code, filename)
        i <- i+1
    }
}

writeLines(ff_in, con = "analysis/hpc/dst_model_input_filenames.txt")
