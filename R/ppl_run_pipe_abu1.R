#' Run abundance indicators pipeline for a site
#'
#' @param site Code for site of interest.
#' @param year Year to run to the pipeline for
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambJAGS}
#' @param steps A character vector containing the steps of the pipeline to run.
#' Can contain: "data", "fit", "summ". Defaults to all of them.
#' @param ... Other arguments passed on to other functions
#'
#' @return
#' @export
#'
#' @examples
ppl_run_pipe_abu1 <- function(site, year, config,
                              steps = c("data", "fit", "summ"), ...){

    # Prepare covariates if necessary
    if("data" %in% steps){

        library(rgee)
        ee_check()
        ee_Initialize(drive = TRUE)

        ppl_data_ssm(site_id, year, config)

        detach("package:rgee", unload = TRUE)

    }

    # Load counts
    counts <- readRDS(file.path(config$data_dir, paste(site, year, "visit_covts.rds", sep = "_")))

    # Loop through selected species
    for(i in seq_along(config$species)){

        sp_code <- config$species[i]

        # Species name
        sp_name <- counts %>%
            dplyr::filter(SppRef == sp_code) %>%
            dplyr::mutate(name = paste(Common_species, Common_group)) %>%
            dplyr::pull(name) %>%
            unique()

        print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

        # Fit model
        if("fit" %in% steps){
            ppl_fit_ssm_model(sp_code, site, year, config, ...)
        }

        # Summary and outputs
        if("summ" %in% steps){
            ppl_summarize_ssm(sp_code, site, year, config, ...)
        }
    }
}
