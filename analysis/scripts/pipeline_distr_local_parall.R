library(BIRDIE)
# library(dplyr)
# library(occuR)

rm(list = ls())

test_years <- c(2012, 2013, 2014, 2015)
run_modules <- 1
ncores <- 2
parall <- ncores != 1
annotate <- FALSE

for(y in seq_along(test_years)){

    year <- test_years[y]
    config <- configPreambOccu(year = year, dur = 3, dim_grid = 10, server = FALSE)

    # Create aux indices for parallel computing
    if(parall){
        seq_spp <- seq_along(config$species)
        keep <- seq(min(seq_spp), max(seq_spp), ncores)

        if((length(seq_spp) %% ncores) == 0){
            ppll <- rep(TRUE, length(keep))
        } else {
            ppll <- c(rep(TRUE, length(keep)-1), FALSE)
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

        # Annotate data if necessary
        if(annotate && i == 1){

            message(paste0("Annotating sites/visits species ", paste(sp_codes[1], collapse = ", "), " (", i, " of ", length(keep), ")"))

            ppl_create_site_visit(sp_code = sp_codes[1],
                                  year = year,
                                  force_gee_dwld = TRUE,
                                  save_occu_data = TRUE,
                                  overwrite_occu_data = c("site", "visit", "det"),
                                  config = config,
                                  force_abap_dwld = TRUE,
                                  monitor = TRUE)
        }

        message(paste0("Working on species ", paste(sp_codes, collapse = ", "), " (", i, " of ", length(keep), ")"))

        if(parall){
            future::plan("multisession", workers = ncores)
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
