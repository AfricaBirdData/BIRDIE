library(BIRDIE)
# library(dplyr)
# library(occuR)

rm(list = ls())

test_years <- c(2012, 2013, 2014, 2015)
run_modules <- 1
ncores <- 4
parall <- ncores != 1
annotate <- FALSE
run_data_prep <- TRUE

for(y in seq_along(test_years)){

    year_sel <- test_years[y]


    # Preamble ----------------------------------------------------------------

    # Create config object
    config <- configPreambOccu(year = year_sel, dur = 3,
                               occ_mod = c("log_dist_coast", "watext", "log_watext", "watrec", "ndvi", "elev",
                                           "prcp", "tdiff", "watext:watrec"),
                               det_mod = c("(1|obs_id)", "(1|site_id)", "log_hours", "prcp", "tdiff", "cwac"),
                               fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "wetext_2018","wetcon_2018",
                                              "dist_coast", "elev"),
                               package = "spOccupancy",
                               server = FALSE)

    # Create log (we do this once for each period defined by config$years)
    if(y == 1){
        createLog(config, log_file = NULL, date_time = NULL, species = NA, model = NA,
                  year = NA, data = NA, fit = NA, diagnose = NA, summary = NA,
                  package = NA, notes = "Log file created")
    }

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


    # Annotate data with GEE --------------------------------------------------

    # Annotate data if necessary (we also only annotate data once each period defined by config$years)
    if(annotate){

        message(paste("Annotating sites/visits for years", config$year_range))

        ppl_run_pipe_dst1(sp_code = config$species[1],
                          sp_name = "sp_name",
                          year = year_sel,
                          config = config,
                          steps = c("data"),
                          force_gee_dwld = TRUE,
                          monitor_gee = FALSE,
                          force_site_visit = TRUE,
                          force_abap_dwld = FALSE,
                          spatial = FALSE,
                          print_fitting = FALSE)
    }


    # RUN PIPELINE ------------------------------------------------------------

    for(k in seq_along(keep)){

        # Index handling for parallel computing
        idx <- keep[k]
        j <- ppll[k]
        sp_codes <- config$species[idx:(idx+j-1)]


        # 1. Run data preparation routines in series ------------------------------

        message(paste0("Preparing data for species ", paste(sp_codes, collapse = ", "), " (", k, " of ", length(keep), ")"))

        # We need to run through all species to prepare detection data, but not site and visit data.
        # That is why force_site_visit = FALSE below, for all species but one
        for(i in seq_along(sp_codes)){

            sp_code <- sp_codes[i]

            # Species name
            sp_name <- BIRDIE::barberspan %>%
                dplyr::filter(SppRef == sp_code) %>%
                dplyr::mutate(name = paste(Common_species, Common_group)) %>%
                dplyr::mutate(name = gsub(" NA|NA ", "", name)) %>% # in case there are NAs in species or group
                dplyr::pull(name) %>%
                unique()

            message(paste("Preparing data for species", sp_code))

            out_dst1 <- ppl_run_pipe_dst1(sp_code = sp_code,
                                          sp_name = sp_name,
                                          year = year_sel,
                                          config = config,
                                          steps = c("data"),
                                          force_gee_dwld = FALSE,
                                          monitor_gee = FALSE,
                                          force_site_visit = ifelse(i == 1, TRUE, FALSE),
                                          force_abap_dwld = FALSE,
                                          spatial = FALSE,
                                          print_fitting = FALSE)

            if(out_dst1 == 1){
                next
            }
        }


        # 2. Run model fitting in parallel ----------------------------------------

        message(paste0("Fitting models for species ", paste(sp_codes, collapse = ", "), " (", k, " of ", length(keep), ")"))

        if(parall){
            future::plan("multisession", workers = ppll[i])
        }

        for(t in seq_along(config$years)){

            year_sel <- config$years[t]

            furrr::future_map(sp_codes, ~pipe_prll_fit(.x, year_sel, .spatial = FALSE, config),
                              .options = furrr::furrr_options(seed = TRUE,
                                                              packages = c("BIRDIE", "spOccupancy")))

        }

        if(parall){
            future::plan("sequential")
        }
    }
}
