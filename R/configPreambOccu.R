#' Prepare server or local occupancy modelling preamble
#'
#' @description Set basic variables to run BIRDIE scripts locally or remotely.
#' @param year Year of interest.
#' @param dur Temporal coverage of the analysis in years. `year` will be the last year
#' covered by the analysis.
#' @param occ_mod A character vector with the names of the variables to include in the occupancy
#' process in the occupancy model. Random effects and interactions are specified as in
#' \code{\link[lme4]{lmer}}. Note that only second order interactions are accepted at the moment
#' (i.e., interactions of two variables).
#' @param det_mod A character vector with the names of the variables to include in the detection
#' process in the occupancy model. Random effects and interactions are specified as in
#' \code{\link[lme4]{lmer}}. Note that only second order interactions are accepted at the moment.
#' (i.e., interactions of two variables).
#' @param package A character string with the name of the package that should be used for
#' fitting occupancy models. Currently: "spOccupancy", "occuR".
#' @param fixed_vars A character vector with the names of the variables included in the
#' occupancy model that don't change over time.
#' @param dim_grid This was for occuR and it is not doing anything at the moment.
#' An integer giving the dimension of the grid used for spatial effects.
#' This dimension gives the number `k` see \code{\link[mgcv]{choose.k}}. Defaults
#' to 10.
#' @param server Logical. If TRUE the preamble is prepared to run remotely,
#' otherwise it is prepared to run locally.
#' @param data_dir Path to data directory. There are a few inputs to the pipeline
#' that it doesn't generate itself, such as some environmental layers that are not
#' on Google Earth Engine. Those would be stored here.
#' @param out_dir Path to output directory. Pipeline outputs will be stored here,
#' including intermediate outputs, so most what we need is here.
#'
#' @return A list of elements to configure the pipeline to process occupancy models
#' @export
#'
#' @examples
#' configPreambOccu(year = 2010, dur = 3, dim_grid = 20, server = TRUE)
configPreambOccu <- function(year, dur, occ_mod, det_mod, package = "spOccupancy",
                             fixed_vars, server, dim_grid = 10, data_dir = NULL,
                             out_dir = NULL){

    if(server){

        # Define data and output directories
        if(is.null(data_dir)){
            data_dir <- "/home/birdie/analysis/data"
        }

        if(is.null(out_dir)){
            out_dir <- "/drv_birdie/birdie_ftp"
        }

        # Define species to fit models to
        species <- unique(BIRDIE::barberspan$SppRef) # For now, we want to select species present at Barberspan

        # Remove partially identified species
        species <- species[species < 10000]

    } else {

        # Define data and output directories
        if(is.null(data_dir)){
            data_dir <- "analysis/data"
        }
        if(is.null(out_dir)){
            out_dir <- "analysis/out_nosync"
        }

        # Define species to fit models to
        species <- c(4, 6, 41, 235, 240)

    }

    # Define a range of years covered by the occupancy model
    year_range <- c(year - dur + 1, year)
    years_ch <- paste(substring(as.character(year_range), 3, 4), collapse = "_")
    years <- year_range[1]:year_range[2]

    list(server=server, data_dir=data_dir, out_dir=out_dir, species=species, dim_grid=dim_grid,
         year=year, dur=dur, year_range=year_range, years_ch=years_ch, years=years,
         occ_mod=occ_mod, det_mod=det_mod, package=package, fixed_vars=fixed_vars)

}
