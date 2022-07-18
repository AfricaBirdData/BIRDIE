#' Set a standard path for a pipeline output file associated with a species
#'
#' @description The datapipeline produces many output files. This function
#' creates an standard output file path for those files that are associated with
#' a particular species. There might be instances where we need name that
#'  deviates from the standard. We need to handle these exceptions
#' case by case. Those files that are not associated with a particular species
#' follow different standards.
#' @param prefix A character string with the prefix that will appear at the
#' beginning of the name. This is what distinguishes files within the species
#' directory.
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambJAGS} and \link{configPreambOccuR}.
#' @param sp_code SAFRING reference number of the species we want to analyze.
#' @param ext The extension of the output file. Note that we must add the
#' trailing '.' to the extension.
#'
#' @return A character string with the path to/for a pipeline output file
#' @export
#'
#' @examples
#' \dontrun{
#' config <- configPreambJAGS(2017, server = FALSE)
#' sp_code <- 4
#'
#' # Set a file path for a state-space model output
#' setSpOutFilePath("ssm_fit", config, sp_code, ".rds)
#'
#' # Set a file path for model-ready data
#' setSpOutFilePath("model_data", config, sp_code, ".csv)
#' }
setSpOutFilePath <- function(prefix, config, sp_code, ext){
    file.path(config$out_dir, sp_code, paste0(prefix, "_", config$years_ch, "_", sp_code, ext))
}
