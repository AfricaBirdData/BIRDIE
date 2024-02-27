#' Configure pipeline parameters
#'
#' @description Set basic variables to run the BIRDIE pipeline locally or remotely.
#' @param year Year of interest.
#' @param dur Temporal coverage of the analysis in years. `year` will be the last year
#' covered by the analysis.
#' @param region A character string with the region we want to run the pipeline for.
#' Currently only "South Africa" and "Kenya" are supported.
#' @param module A character string defining the module the pipeline should run.
#' At the moment it can be one of `c("dst", "abu")` for distributions and abundance
#' respectively.
#' @param occ_mod A character vector with the names of the variables to include in the occupancy
#' process in the occupancy model. Random effects and interactions are specified as in
#' \code{\link[lme4]{lmer}}. Note that only second order interactions are accepted at the moment
#' (i.e., interactions of two variables).
#' @param det_mod A character vector with the names of the variables to include in the detection
#' process in the occupancy model. Random effects and interactions are specified as in
#' \code{\link[lme4]{lmer}}. Note that only second order interactions are accepted at the moment.
#' (i.e., interactions of two variables).
#' @param fixed_vars A character vector with the names of the variables included in the
#' occupancy model that don't change over time.
#' @param mod_file Name of file containing model, with out path to directory.
#' Directory is specified in `mod_dir`. This is typically used for JAGS or Stan
#' where models are written on an external file.
#' @param server Logical. If TRUE the preamble is prepared to run remotely,
#' otherwise it is prepared to run locally.
#' @param data_dir Path to data directory. There are a few inputs to the pipeline
#' that it doesn't generate itself, such as some environmental layers that are not
#' on Google Earth Engine. Those would be stored here.
#' @param mod_dir Path to directory where models are saved.
#' @param out_dir Path to output directory. Pipeline outputs will be stored here,
#' including intermediate outputs, so most what we need is here.
#' @param package A character string with the name of the package that should be used for
#' fitting occupancy or state-space models.
#'
#' @return A list of parameters that will be passed on to other functions in the
#' pipeline.
#' @export
#'
#' @examples
#' config <- configPipeline(
#'     year = 2021,
#'     dur = 29,
#'     mod_file = "cwac_ssm_two_season_mean_rev.R",
#'     package = "jagsUI",
#'     data_dir = NULL,
#'     out_dir = NULL,
#'     server = FALSE
#'     )
configPipeline <- function(year, dur, region = c("southafrica", "kenya"), module = c("dst", "abu"),
                           occ_mod = NULL, det_mod = NULL, fixed_vars = NULL,
                           mod_file = NULL, server = FALSE, data_dir = NULL,
                           mod_dir = NULL, out_dir = NULL, package = NULL){

    region <- gsub(" ", "", region)
    region <- tolower(region)
    region <- match.arg(region)

    region_short <- dplyr::case_when(region == "southafrica" ~ "ZA",
                                     region == "kenya" ~ "KE")

    module <- match.arg(module)

    # Check interaction terms in model specification. The order in the interactions
    # variables must be the same in the interaction term and in the variable list
    if(!is.null(occ_mod)){
        occ_intrcs <- occ_mod[grepl("\\:", occ_mod)]
        occ_intrcs_vars <- strsplit(occ_intrcs, ":")
        for(i in seq_along(occ_intrcs_vars)){
            occ_mod[occ_mod %in% occ_intrcs_vars[[i]]] <- occ_intrcs_vars[[i]]
        }
    }

    if(!is.null(det_mod)){
        det_intrcs <- det_mod[grepl("\\:", det_mod)]
        det_intrcs_vars <- strsplit(det_intrcs, ":")
        for(i in seq_along(det_intrcs_vars)){
            det_mod[det_mod %in% det_intrcs_vars[[i]]] <- det_intrcs_vars[[i]]
        }
    }


    # Define species to fit models to
    species <- sort(unique(BIRDIE::waterbirds$SppRef)) # For now, we want to select species present at Barberspan

    # Remove species without a code
    species <- species[!is.na(species)]


    if(server){

        # Define data and output directories
        if(is.null(data_dir)){
            data_dir <- "/home/birdie/analysis/data"
        }

        if(is.null(out_dir)){
            out_dir <- "/drv_birdie/birdie_ftp"
        }

        if(is.null(mod_dir)){
            mod_dir <- "/drv_birdie/Working/git/BIRDIE/analysis/models"
        }

    } else {

        # Define data and output directories
        if(is.null(data_dir)){
            data_dir <- "analysis/data"
        }

        if(is.null(out_dir)){
            out_dir <- "analysis/output"
        }

        if(is.null(mod_dir)){
            mod_dir <- "analysis/models"
        }

        # Define species to fit models to
        # species <- c(4, 6, 41, 235, 240)

    }

    # Define a range of years covered by the model
    # Define a range of years covered by the occupancy model
    year_range <- c(year - dur + 1, year)
    years_ch <- paste(substring(as.character(year_range), 3, 4), collapse = "_")
    years <- year_range[1]:year_range[2]

    list(server=server, region=region_short, module=module, data_dir=data_dir, out_dir=out_dir, mod_dir=mod_dir, mod_file=mod_file,
         species=species, year=year, dur=dur, year_range=year_range, years_ch=years_ch, years=years,
         occ_mod=occ_mod, det_mod=det_mod, package=package, fixed_vars=fixed_vars)

}
