library(BIRDIE)

rm(list = ls())


# Configuration -----------------------------------------------------------

# We are currently working with several detection models
det_mods <- list(det_mod1 = c("(1|obs_id)", "log_hours", "prcp", "tdiff", "cwac"),
                 det_mod2 = c("(1|obs_id)", "(1|site_id)", "log_hours", "prcp", "tdiff", "cwac"))

# Configure pipeline
config <- configPipeline(year = 2010,
                         dur = 3,
                         module = "dst",
                         occ_mod = c("log_dist_coast", "elev", "log_hum.km2", "wetcon",
                                     "watrec", "watext", "log_watext", "watrec:watext",
                                     "ndvi", "prcp", "tdiff"),
                         det_mod = det_mods$det_mod1,
                         fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "wetext_2018","wetcon_2018",
                                        "dist_coast", "elev"),
                         package = "spOccupancy",
                         data_dir = NULL,     # this might have to be adapted?
                         out_dir = NULL,     # this might have to be adapted?
                         server = FALSE)


for(i in seq_along(config$species)){

    # Select one species code
    sp_code <- config$species[i]


    # Create logs -------------------------------------------------------------

    # Create log?
    if(i == 1){
        createLog(config, log_file = NULL, date_time = NULL, species = NA, model = NA,
                  year = NA, data = NA, fit = NA, diagnose = NA, summary = NA,
                  package = NA, notes = "Log file created")
    }

    message(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))


    # Run modules -------------------------------------------------------------

    for(t in seq_along(config$years)){

        year_sel <- config$years[t]

        out_dst1 <- ppl_run_pipe_dst1(sp_code = sp_code,
                                      year = year_sel,
                                      config = config,
                                      steps = c("data", "fit", "diagnose", "summary"),
                                      force_gee_dwld = FALSE,
                                      monitor_gee = TRUE,
                                      force_site_visit = TRUE,
                                      force_abap_dwld = FALSE,
                                      spatial = FALSE,
                                      print_fitting = TRUE)

        message(paste("Pipeline DST1 status =", out_dst1))

        if(out_dst1 != 0){
            next
        }
    }
}
