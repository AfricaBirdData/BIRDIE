library(BIRDIE)
# library(dplyr)
# library(occuR)

rm(list = ls())

test_years <- c(2012, 2013, 2014, 2015)
run_modules <- 1
ncores <- 4
parall <- ncores != 1
annotate <- FALSE

for(y in seq_along(test_years)){

    year <- test_years[y]
    config <- configPreambOccu(year = year, dur = 3,
                               occ_mod = c("log_dist_coast", "watext", "log_watext", "watrec", "ndvi", "elev",
                                           "prcp", "tdiff", "watext:watrec"),
                               det_mod = c("(1|obs_id)", "(1|site_id)", "log_hours", "prcp", "tdiff", "cwac"),
                               fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "log_dist_coast", "elev"),
                               package = "spOccupancy",
                               server = TRUE)

    createLog(config, logfile = NULL, date_time = NULL, species = NA, model = NA,
              year = NA, data = NA, fit = NA, diagnose = NA, summary = NA,
              package = NA, notes = "Log file created")

    # Create aux indices for parallel computing
    if(parall){
        seq_spp <- seq_along(config$species)
        keep <- seq(min(seq_spp), max(seq_spp), ncores)

        last_run <- (length(seq_spp) %% ncores)

        if(last_run == 1){
            ppll <- c(rep(ncores, length(keep)-1), 1)
        } else if(last_run == 0){
            ppll <- rep(ncores, length(keep))
        } else {
            ppll <- c(rep(ncores, length(keep)-1), last_run)
        }
    } else {
        keep <- seq_along(config$species)
        ppll <- rep(1, length(keep))
    }

    for(i in seq_along(keep)){

        # Index handling for parallel computing
        idx <- keep[i]
        j <- ppll[i]
        sp_codes <- config$species[idx:(idx+j-1)]

        # Annotate data if necessary
        if(annotate && i == 1){

            message(paste0("Annotating sites/visits species ", paste(sp_codes[1], collapse = ", "), " (", i, " of ", length(keep), ")"))

            ppl_create_site_visit(sp_code = sp_codes[1],
                                  force_gee_dwld = TRUE,
                                  save_occu_data = TRUE,
                                  overwrite_occu_data = c("site", "visit", "det"),
                                  config = config,
                                  force_abap_dwld = TRUE,
                                  monitor = TRUE)
        }

        message(paste0("Working on species ", paste(sp_codes, collapse = ", "), " (", i, " of ", length(keep), ")"))

        if(parall){
            future::plan("multisession", workers = ppll[i])
        }

        for(t in seq_along(config$years)){

            year_sel <- config$years[t]

            furrr::future_map(sp_codes, ~pipe_prll_fit(.x, year_sel, .spatial = FALSE, config))

        }

        if(parall){
            future::plan("sequential")
        }

    }
}
