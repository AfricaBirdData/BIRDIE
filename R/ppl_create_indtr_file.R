#' Create indicator storage files
#'
#' @param sp_code SAFRING code of the species to run the pipeline for
#' @param year Year to run to the pipeline for
#' @param overwrite_indtr Logical. If TRUE, existing files in directories
#' corresponding to the species in config$species will be overwritten.
#' @param update_file A character string with the file path to a previous indicator
#' file that should serve as a basis for the new indicator file. Defaults to NULL,
#' and so a clean new file would be created.
#'
#' @return
#' @export
#'
#' @examples
ppl_create_indtr_file <- function(sp_code, year, overwrite_indtr, update_file = NULL){

    # Create indicator file path
    indtr_file <- file.path(config$out_dir, sp_code, paste0("indtr_dst_", sp_code, "_", year, ".csv"))

    # Check if file exists
    if(!file.exists(indtr_file) | (file.exists(indtr_file) & overwrite_indtr)){

        if(is.null(update_file)){

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

        } else {

            indtr <- read.csv(update_file)

        }

        # Save
        write.csv(indtr, indtr_file, row.names = FALSE)

    }
}
