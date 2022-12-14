#' @export
#' @example
#' \dontrun{}
#'
#'
#' Fit occupancy models in parallel
#'
#' @description Wrapper around ppl_run_pipe_dst1 (fit step) for furrr parallel calls
#' @param .sp_code Species SAFRING code
#' @param .year Year to run to the pipeline for
#' @param .spatial Whether a spatial model should be fit. Defaults to FALSE.
#' @param .config Config object from \link{configPreambOccu}
#'
#' @keywords internal
#' @noRd
pipe_prll_fit <- function(.sp_code, .year, .spatial = FALSE, .config){

    # Species name
    sp_name <- BIRDIE::barberspan %>%
        dplyr::filter(SppRef == .sp_code) %>%
        dplyr::mutate(name = paste(Common_species, Common_group)) %>%
        dplyr::mutate(name = gsub(" NA|NA ", "", name)) %>% # in case there are NAs in species or group
        dplyr::pull(name) %>%
        unique()

    # Pipeline module 1
    out_dst1 <- ppl_run_pipe_dst1(sp_code = .sp_code,
                                  sp_name = sp_name,
                                  year = .year,
                                  config = .config,
                                  steps = c("data", "fit"),
                                  force_gee_dwld = FALSE,
                                  force_abap_dwld = FALSE,
                                  save_occu_data = TRUE,
                                  overwrite_occu_data = c("det"),
                                  scale_vars_occur = list(visit = NULL,
                                                          site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")),
                                  spatial = .spatial,
                                  print_fitting = FALSE,
                                  verbose = TRUE,
                                  monitor = FALSE)

    message(paste("Pipeline DST1 status =", out_dst1))

}
