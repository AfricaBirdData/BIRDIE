library(BIRDIE)

rm(list = ls())

# We are currently working with several detection models
det_mods <- list(det_mod1 = c("(1|obs_id)", "log_hours", "prcp", "tdiff", "cwac"),
                 det_mod2 = c("(1|obs_id)", "(1|site_id)", "log_hours", "prcp", "tdiff", "cwac"))

# Configure pipeline
ncores <- 19
parall <- ncores != 1
annotate <- FALSE
prep_site_visit <- FALSE

config <- configPipeline(year = 2019,
                         dur = 3,
                         module = "dst",
                         occ_mod = c("log_dist_coast", "elev", "log_hum.km2", "wetcon",
                                     "watrec", "watext", "log_watext", "watrec:watext",
                                     "ndvi", "prcp", "tdiff"),
                         det_mod = det_mods$det_mod2,
                         fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "wetext_2018","wetcon_2018",
                                        "dist_coast", "elev"),
                         package = "spOccupancy",
                         data_dir = "/home/crvfra001/birdie/data",
                         out_dir = "/scratch/crvfra001/birdie/output",
                         server = TRUE)

# Create log (we do this once for each period defined by config$years)
createLog(config, log_file = NULL, date_time = NULL, species = NA, model = NA,
          year = NA, data = NA, fit = NA, diagnose = NA, summary = NA,
          package = NA, notes = "Log file created")


# Annotate data with GEE --------------------------------------------------

# Annotate data if necessary (we also only annotate data once each period defined by config$years)
if(annotate){

    message(paste("Annotating sites/visits for years", config$year_range))

    ppl_run_pipe_dst1(sp_code = config$species[1],
                      sp_name = "sp_name",
                      year = config$years[1],  # This year could be anything really
                      config = config,
                      steps = c("data"),
                      force_gee_dwld = TRUE,
                      monitor_gee = TRUE,
                      force_site_visit = FALSE,
                      force_abap_dwld = FALSE,
                      spatial = FALSE,
                      print_fitting = FALSE)
}


# RUN PIPELINE ------------------------------------------------------------

# Add some extra species for Alan
# sp_codes <- c(config$species, 566, 463, 220, 4127, 461, 764, 653, 626, 613, 4125, 1037, 486, 479)
sp_codes <- c(95, 91, 96, 89, 88, 102, 90, 94, 98, 99, 97, 100, 274, 233, 237, 238, 235, 245, 246, 281, 231, 228, 288, 289, 287, 290, 291, 296, 298, 305, 304, 269, 270, 256, 258, 264, 250, 253, 263, 268, 52, 61, 58, 59, 60, 64, 55, 56, 63, 54, 57, 62, 74, 73, 75, 77, 76, 42, 41, 48, 50, 47, 86, 87, 5, 4, 6, 72, 81, 83, 84, 85, 397, 394, 395, 149, 167, 216, 214, 203, 210, 197, 208, 212, 764, 4125)
# sp_codes <- sp_codes[1:100]

# 1. Run data preparation routines in series ------------------------------

# We need to run through all species to prepare detection data, but not site and visit data.
# That is why force_site_visit = FALSE below, for all species but one
if(prep_site_visit){

    for(i in seq_along(sp_codes)){

        sp_code <- sp_codes[i]

        message(paste("Preparing data for species", sp_code))

        out_dst1 <- ppl_run_pipe_dst1(sp_code = sp_code,
                                      year = config$years[1],  # This year could be anything really
                                      config = config,
                                      steps = c("data"),
                                      force_gee_dwld = FALSE,
                                      monitor_gee = FALSE,
                                      force_site_visit = ifelse(i == 1, TRUE, FALSE),
                                      force_abap_dwld = FALSE,
                                      spatial = FALSE,
                                      print_fitting = FALSE)

        if(out_dst1 != 0){
            next
        }
    }
}


# 2. Run model fitting in parallel ----------------------------------------

message(paste0("Processing species ", paste(sp_codes, collapse = ", ")))

if(parall){
    future::plan("multisession", workers = ncores)
}

for(t in seq_along(config$years)){

    if(t %in% c(1, 2)) next

    year_sel <- config$years[t]

    furrr::future_map(sp_codes, ~pipe_prll_fit(.x, year_sel, .spatial = FALSE, config,
                                               .steps = c("fit", "diagnose"), time_limit = 12*3600),
                      .options = furrr::furrr_options(seed = TRUE,
                                                      packages = c("BIRDIE", "spOccupancy")))

}

if(parall){
    future::plan("sequential")
}
