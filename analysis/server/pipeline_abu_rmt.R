library(BIRDIE)
library(dplyr)

rm(list = ls())

year <- 2006

config <- configPreambJAGS(year = year, server = TRUE)

site <- 26352535

for(i in seq_along(config$species)){

    sp_code = config$species[i]

    # Species name
    sp_name <- BIRDIE::barberspan %>%
        dplyr::filter(SppRef == sp_code) %>%
        mutate(name = paste(Common_species, Common_group)) %>%
        pull(name) %>%
        unique()

    print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

    ppl_run_pipe_abu(sp_code = sp_code,
                     site = site,
                     year = year,
                     config = config,
                     steps = c("fit", "summ"))

}
