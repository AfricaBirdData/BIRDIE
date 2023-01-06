library(BIRDIE)
# library(dplyr)
# library(occuR)

rm(list = ls())

test_years <- c(2012, 2013, 2014, 2015)

run_modules <- 1

y=1; i=1; t=1

for(y in seq_along(test_years)){

    year <- test_years[y]
    config <- configPreambOccu(year = year, dur = 3,
                               occ_mod = c("log_dist_coast", "watext", "log_watext", "watrec", "ndvi", "elev",
                                            "prcp", "tdiff", "watext:watrec"),
                               det_mod = c("(1|obs_id)", "(1|site_id)", "log_hours", "prcp", "tdiff", "cwac"),
                               fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "dist_coast", "elev"),
                               server = FALSE)

    createLog(config, logfile = NULL, date_time = NULL, species = NA, model = NA, data = NA, fit = NA,
              diagnose = NA, summary = NA, notes = "Log file created")

    for(i in seq_along(config$species)){

        sp_code <- config$species[i]

        # Species name
        sp_name <- BIRDIE::barberspan %>%
            dplyr::filter(SppRef == sp_code) %>%
            dplyr::mutate(name = paste(Common_species, Common_group)) %>%
            dplyr::mutate(name = gsub(" NA|NA ", "", name)) %>% # in case there are NAs in species or group
            dplyr::pull(name) %>%
            unique()

        message(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

        if(1 %in% run_modules){

            for(t in seq_along(config$years)){

                year_sel <- config$years[t]

                out_dst1 <- ppl_run_pipe_dst1(sp_code = sp_code,
                                              sp_name = sp_name,
                                              year = year_sel,
                                              config = config,
                                              steps = c("summary"),
                                              force_gee_dwld = FALSE,
                                              force_abap_dwld = FALSE,
                                              save_occu_data = TRUE,
                                              overwrite_occu_data = c("site", "visit", "det"),
                                              scale_vars_occur = list(visit = NULL,
                                                                      site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")),
                                              spatial = FALSE,
                                              print_fitting = FALSE,
                                              verbose = TRUE,
                                              monitor = TRUE)

                message(paste("Pipeline DST1 status =", out_dst1))

                if(out_dst1 == 1){
                    next
                }

            }

        }


        if(2 %in% run_modules){

            out_dst2 <- ppl_run_pipe_dst2(sp_code = sp_code,
                                          config = config,
                                          indtr = c("aoo", "daoo"),
                                          create = TRUE,
                                          overwrite_indtr = TRUE,
                                          verbose = TRUE,
                                          scale_vars_occur = list(visit = NULL,
                                                                  site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")))

            message(paste("Pipeline DST2 status =", out_dst2))

            if(out_dst2 == 1){
                next
            }

        }
    }
}
