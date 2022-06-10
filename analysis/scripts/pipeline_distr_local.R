library(BIRDIE)
# library(dplyr)
# library(occuR)

rm(list = ls())

test_years <- c(2012, 2013, 2014)

run_branches <- c(1, 2)


for(y in seq_along(test_years)){
y=2
    year <- test_years[y]
    config <- configPreambOccuR(year = year, dur = 3, dim_grid = 10, server = FALSE)

    #for(i in seq_along(config$species)){
i=1
        sp_code <- config$species[i]

        # Species name
        sp_name <- BIRDIE::barberspan %>%
            dplyr::filter(SppRef == sp_code) %>%
            dplyr::mutate(name = paste(Common_species, Common_group)) %>%
            dplyr::mutate(name = gsub(" NA|NA ", "", name)) %>% # in case there are NAs in species or group
            dplyr::pull(name) %>%
            unique()

        print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

        if(1 %in% run_branches){

            out_dst1 <- ppl_run_pipe_dst1(sp_code = sp_code,
                                          sp_name = sp_name,
                                          year = year,
                                          config = config,
                                          steps = c("data", "fit", "summ"),
                                          force_gee_dwld = FALSE,
                                          force_abap_dwld = FALSE,
                                          save_occu_data = TRUE,
                                          overwrite_occu_data = c("site", "visit", "det"),
                                          scale_vars_occur = list(visit = NULL,
                                                                  site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")),
                                          print_fitting = FALSE,
                                          verbose = TRUE)

            message(paste("Pipeline DST1 status =", out_dst1))

            if(out_dst1 == 1){
                next
            }

        }


        if(2 %in% run_branches){

            out_dst2 <- ppl_run_pipe_dst2(sp_code = sp_code,
                                          config = config,
                                          indtr = c("aoo", "daoo"),
                                          create = FALSE,
                                          overwrite_indtr = TRUE,
                                          verbose = TRUE,
                                          scale_vars_occur = list(visit = NULL,
                                                                  site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")))

            message(paste("Pipeline DST2 status =", out_dst2))

            if(out_dst2 == 1){
                next
            }

        }
        #}
}
