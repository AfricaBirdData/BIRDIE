#' Create indicator storage files
#'
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambOccuR}
#' @param overwrite_indtr Logical. If TRUE, existing files in directories
#' corresponding to the species in config$species will be overwritten.
#'
#' @return
#' @export
#'
#' @examples
ppl_create_indtr_files <- function(config, overwrite_indtr){

    for(i in seq_along(config$species)){

        # Select species
        sp_code = config$species[i]

        # Create indicator file path
        indtr_path <- file.path(config$fit_dir, sp_code, paste0("indtr_", sp_code, ".csv"))

        # Check if file exists
        if(!file.exists(indtr_path) | (file.exists(indtr_path) & overwrite_indtr)){

            # Create empty data frame
            indtr <- data.frame(species = character(),
                                indicator = character(),
                                start_date = character(),
                                end_date = character(),
                                term = character(),
                                estimate = numeric(),
                                st_dev = numeric(),
                                lb95 = numeric(),
                                ub95 = numeric(),
                                opt = numeric())

            # Save
            write.csv(indtr, indtr_path, row.names = FALSE)

        }
    }
}
