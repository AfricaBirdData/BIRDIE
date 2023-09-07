#' Select species based on model diagnostics
#'
#' @description The datapipeline produces one occupancy diagnostic file for each
#' species and year. This function combines all these files and returns those
#' species with convergence or goodnes-of-fit issues.
#' @param config A list with pipeline configuration parameters.
#' See \link{configPipeline}.
#' @param year The year for which diagnostics are required.
#' @param sp_codes SAFRING reference numbers of the species we want diagnostics for.
#' @param module Either "dst" if we want occupancy diagnostics or "abu" if we want
#' state-space model diagnostics.
#'
#' @return A list with two elements: $no_converge, contains codes of all species
#' with convergence issues, $bad_fit, contains codes of species with goodness-of-fit
#' issues (small Bayesian p-value).
#' @export
#'
#' @examples
#' \dontrun{
#' config <- configPipeline(year = 2010,
#'                          dur = 3,
#'                          occ_mod = c("log_dist_coast", "elev", "log_hum.km2", "wetcon",
#'                                      "watrec", "watext", "log_watext", "watext:watrec",
#'                                      "ndvi", "prcp", "tdiff"),
#'                          det_mod = c("(1|obs_id)", "log_hours", "prcp", "tdiff", "cwac"),
#'                          fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "wetext_2018","wetcon_2018",
#'                                         "dist_coast", "elev"),
#'                          package = "spOccupancy",
#'                          data_dir = "analysis/hpc/imports",
#'                          out_dir = "analysis/hpc/imports",
#'                          server = TRUE)
#' sp_codes <- config$species
#'
#' selectSppFromDiag(config, sp_codes, 2008)
#' }
selectSppFromDiag <- function(config, sp_codes, year, module = c("dst", "abu")){

    module <- match.arg(module)

    # Combine species diagnostics
    if(module == "dst"){

        diags <- combineOccuDiags(config, sp_codes, year)

        # Create output
        out <- vector("list")

        # Which species have non-convergent parameters
        out$no_converge <- diags %>%
            dplyr::filter(nc_pars != 0) %>%
            dplyr::pull(sp)

        # Which species have significant Bayesian p-value
        out$bad_fit <- diags %>%
            dplyr::filter(bayes_p < 0.05) %>%
            dplyr::pull(sp)

    } else if(module == "abu"){

        diags <- combineAbuDiags(config, sp_codes, year)

        # Create output
        out <- vector("list")

        # Which species have non-convergent parameters
        out$no_converge <- diags %>%
            dplyr::filter(nc_pars != 0) %>%
            dplyr::pull(sp)

        # Which species have significant Bayesian p-value
        out$bad_fit <- diags %>%
            dplyr::filter(Tmean < 0.05, Tmean > 0.95,
                          Tsd < 0.05, Tsd > 0.95,
                          abs(Tdiff) > 0.5) %>%
            dplyr::pull(sp)

    }

    return(out)

}
