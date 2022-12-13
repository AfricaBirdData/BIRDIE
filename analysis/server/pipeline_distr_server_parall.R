library(BIRDIE)
# library(dplyr)
# library(occuR)

rm(list = ls())

test_years <- c(2012, 2013, 2014, 2015)
run_modules <- 1
ncores <- 2
parall <- ncores != 1

for(y in seq_along(test_years)){

    year <- test_years[y]
    config <- configPreambOccu(year = year, dur = 3, dim_grid = 10, server = TRUE)

    # Create aux variables for parallel computing
    if(parall){
        seq_spp <- seq_along(config$species)[-1]
        keep <- c(1, seq(min(seq_spp), max(seq_spp), ncores))

        if((length(seq_spp) %% ncores) == 0){
            ppll <- c(FALSE, rep(TRUE, length(keep)-1))
        } else {
            ppll <- c(FALSE, rep(TRUE, length(keep)-2), FALSE)
        }
    } else {
        keep <- seq_along(config$species)
        ppll <- rep(FALSE, length(keep))
    }

    for(i in seq_along(keep)){

        # Index handling for parallel computing
        idx <- keep[i]
        j <- ppll[i]
        sp_codes <- config$species[idx:(idx+j)]

        message(paste0("Working on species ", paste(sp_codes, collapse = ", "), " (", i, " of ", length(keep), ")"))

        pipe_fun <- function(sp_code, .run_modules=run_modules, .config=config, .i=i){

            # Species name
            sp_name <- BIRDIE::barberspan %>%
                dplyr::filter(SppRef == sp_code) %>%
                dplyr::mutate(name = paste(Common_species, Common_group)) %>%
                dplyr::mutate(name = gsub(" NA|NA ", "", name)) %>% # in case there are NAs in species or group
                dplyr::pull(name) %>%
                unique()

            # Pipeline module 1
            if(1 %in% .run_modules){

                out_dst1 <- ppl_run_pipe_dst1(sp_code = sp_code,
                                              sp_name = sp_name,
                                              year = year,
                                              config = .config,
                                              steps = c("data", "fit"),
                                              force_gee_dwld = FALSE,
                                              force_abap_dwld = FALSE,
                                              save_occu_data = TRUE,
                                              overwrite_occu_data = c("site", "visit", "det"),
                                              scale_vars_occur = list(visit = NULL,
                                                                      site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")),
                                              print_fitting = FALSE,
                                              verbose = TRUE,
                                              monitor = TRUE)

                message(paste("Pipeline DST1 status =", out_dst1))

                if(out_dst1 == 1){
                    next
                }

            }

            # Pipeline module 2
            if(2 %in% .run_modules){

                out_dst2 <- ppl_run_pipe_dst2(sp_code = sp_code,
                                              config = .config,
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

        if(parall){
            future::plan("multisession", workers = ncores)
        }

        furrr::future_map(sp_codes, ~pipe_fun(.x))

        if(parall){
            future::plan("sequential")
        }

    }
}
