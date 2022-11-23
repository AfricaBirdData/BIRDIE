#' Prepare server or local occupancy modelling preamble
#'
#' @description Set basic variables to run BIRDIE scripts locally or remotely.
#' @param year Year of interest.
#' @param dur Temporal coverage of the analysis in years. `year` will be the last year
#' covered by the analysis.
#' @param dim_grid This was for occuR and it is not doing anything at the moment.
#' An integer giving the dimension of the grid used for spatial effects.
#' This dimension gives the number `k` see \code{\link[mgcv]{choose.k}}.
#' @param server Logical. If TRUE the preamble is prepared to run remotely,
#' otherwise it is prepared to run locally.
#'
#' @return A list with:
#'   - data_dir: directory where data is retrieved from.
#'   - fit_dir: directory where fitted model and other output is saved to.
#'   - dyear: number of years before and after the year of interest.
#'   - sptemp: spatio-temporal effect of the occupancy model.
#'   - species: species models will be fitted to.
#' @export
#'
#' @examples
#' configPreambOccu(year = 2010, dur = 3, dim_grid = 20, server = TRUE)
configPreambOccu <- function(year, dur, dim_grid, server){

    if(server){

        # Define data and output directories
        data_dir <- "/home/birdie/analysis/data"
        out_dir <- "/drv_birdie/birdie_ftp"

        # Define species to fit models to
        species <- unique(BIRDIE::barberspan$SppRef) # For now, we want to select species present at Barberspan

        # Remove partially identified species
        species <- species[species < 10000]

    } else {

        # Define data and output directories
        data_dir <- "analysis/data"
        out_dir <- "analysis/out_nosync"

        # Define species to fit models to
        species <- c(4, 6, 41, 235, 240)

    }

    # Define a range of years covered by the occupancy model
    year_range <- c(year - dur + 1, year)
    years_ch <- paste(substring(as.character(year_range), 3, 4), collapse = "_")
    years <- year_range[1]:year_range[2]

    list(server=server, data_dir=data_dir, out_dir=out_dir, species=species, dim_grid=dim_grid,
         year=year, dur=dur, year_range=year_range, years_ch=years_ch, years=years)

}
