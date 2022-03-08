#' Prepare server or local preamble
#'
#' @description Set basic variables to run BIRDIE scripts locally or remotely.
#' @param year Year of interest.
#' @param site Code for site of interest.
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
#' configPreambJAGS(year = 2010, server = TRUE)
configPreambJAGS <- function(year, site, server){

    if(server){

        # Define data and output directories
        data_dir <- "/home/birdie/analysis/data"
        mod_file <- "/drv_birdie/Working/git/BIRDIE/analysis/models/cwac_ssm_lat_season.jags"
        data_outdir <- "/drv_birdie/birdie_ftp"
        plot_outdir <- "/drv_birdie/birdie_ftp"

        # Define years to fit
        dyear <- 13

        # Define species to fit models to
        species <- unique(BIRDIE::barberspan$SppRef) # For now, we want to select species present at Barberspan

    } else {

        # Define data and output directories
        data_dir <- "analysis/data"
        mod_file <- "analysis/models/cwac_ssm_lat_season.jags"
        data_outdir <- "analysis/out_nosync"
        plot_outdir <- "analysis/out_nosync"

        # Define years to fit
        dyear <- 24

        # Define species to fit models to
        species <- c(4, 6, 41, 235, 240)

    }

    # Define a range of years covered by the occupancy model
    year_range <- c(year - dyear, year)
    years_ch <- paste(substring(as.character(year_range), 3, 4), collapse = "_")
    years <- year_range[1]:year_range[2]

    list(server=server, data_dir=data_dir, mod_file=mod_file,
         data_outdir=data_outdir, plot_outdir=plot_outdir, species=species,
         year=year, dyear=dyear, year_range=year_range, years_ch=years_ch, years=years)

}
