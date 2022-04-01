#' Create indicator storage files
#'
#' @param sp_code SAFRING code of the species to run the pipeline for
#' @param overwrite_indtr Logical. If TRUE, existing files in directories
#' corresponding to the species in config$species will be overwritten.
#'
#' @return
#' @export
#'
#' @examples
ppl_create_indtr_file <- function(sp_code, overwrite_indtr){

    # Create indicator file path
    indtr_file <- file.path(config$out_dir, sp_code, paste0("indtr_dst_", sp_code, ".csv"))

    # Check if file exists
    if(!file.exists(indtr_file) | (file.exists(indtr_file) & overwrite_indtr)){

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
        write.csv(indtr, indtr_file, row.names = FALSE)

    }
}
