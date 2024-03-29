library(BIRDIE)

rm(list = ls())

# We are currently working with several detection models
det_mods <- list(det_mod1 = c("(1|obs_id)", "log_hours", "prcp", "tdiff", "cwac"),
                 det_mod2 = c("(1|obs_id)", "(1|site_id)", "log_hours", "prcp", "tdiff", "cwac"))

# Configure pipeline
ncores <- 10
parall <- ncores != 1
annotate <- FALSE
prep_site_visit <- FALSE

config <- configPipeline(year = 2020,
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

sp_codes <- config$species

# Add some extra species if necessary
# sp_codes <- c(config$species, 566, 463, 220, 4127, 461, 764, 653, 626, 613, 4125, 1037, 486, 479)

# Run only for some species
# sp_codes <- c(95, 103, 96, 89, 88, 102, 90, 94, 98, 99, 97, 104, 100, 237, 238, 235, 245)

# Split species into groups
# sp_codes <- split(sp_codes, cut(seq_along(sp_codes), 3, labels = FALSE))[[1]]

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

    # Run only for certain year
    # if(t %in% c(2, 3)) next

    year_sel <- config$years[t]

    furrr::future_map(sp_codes, ~pipe_prll_fit(.x, year_sel, .spatial = FALSE, config,
                                               .steps = c("fit", "diagnose"), time_limit = 14*3600),
                      .options = furrr::furrr_options(seed = TRUE,
                                                      packages = c("BIRDIE", "spOccupancy")))


    # If fitting detection model1,
    if(identical(config$det_mod, det_mods$det_mod1)){

        #identify species for which diagnostics failed
        spp <- selectSppFromDiag(config, sp_codes, year_sel)
        bad_spp <- spp$bad_fit

        # Run model with detection model 2 if model 1 diagnostics were not satisfactory
        config$det_mod <- det_mods$det_mod2

        furrr::future_map(bad_spp, ~pipe_prll_fit(.x, year_sel, .spatial = FALSE, config,
                                                  .steps = c("fit", "diagnose"), time_limit = 14*3600),
                          .options = furrr::furrr_options(seed = TRUE,
                                                          packages = c("BIRDIE", "spOccupancy")))
    }

}





if(parall){
    future::plan("sequential")
}
