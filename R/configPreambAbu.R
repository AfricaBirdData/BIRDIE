#' Prepare server or local preamble
#'
#' @description Set basic variables to run BIRDIE scripts locally or remotely.
#' @param year Year of interest.
#' @param server Logical. If TRUE the preamble is prepared to run remotely,
#' otherwise it is prepared to run locally.
#' @param mod_file Name of file containing model, with out path to directory.
#' Directory is specified in `mod_dir`.
#' @param data_dir Path to data directory. There are a few inputs to the pipeline
#' that it doesn't generate itself, such as some environmental layers that are not
#' on Google Earth Engine. Those would be stored here.
#' @param mod_dir Path to directory where models are saved.
#' @param out_dir Path to output directory. Pipeline outputs will be stored here,
#' including intermediate outputs, so most what we need is here.
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
#' configPreambAbu(year = 2010, server = TRUE)
configPreambAbu <- function(year, server, mod_file, data_dir = NULL,
                            mod_dir = NULL, out_dir = NULL){

    if(server){

        # Define data and output directories
        if(is.null(data_dir)){
            data_dir <- "/home/birdie/analysis/data"
        }

        if(is.null(out_dir)){
            out_dir <- "/drv_birdie/birdie_ftp"
        }

        if(is.null(out_dir)){
            mod_dir <- "/drv_birdie/Working/git/BIRDIE/analysis/models"
        }

        # Define years to fit
        dyear <- 28

        # Define species to fit models to
        species <- unique(BIRDIE::barberspan$SppRef) # For now, we want to select species present at Barberspan

    } else {

        # Define data and output directories
        if(is.null(data_dir)){
            data_dir <- "analysis/data"
        }

        if(is.null(out_dir)){
            out_dir <- "analysis/out_nosync"
        }

        if(is.null(mod_dir)){
            mod_dir <- "analysis/models"
        }

        # Define years to fit
        dyear <- 28

        # Define species to fit models to
        species <- c(4, 6, 41, 235, 240)

    }

    # Define a range of years covered by the occupancy model
    year_range <- c(year - dyear, year)
    years_ch <- paste(substring(as.character(year_range), 3, 4), collapse = "_")
    years <- year_range[1]:year_range[2]

    list(server=server, data_dir=data_dir, out_dir=out_dir, mod_dir=mod_dir, mod_file=mod_file,
         species=species, year=year, dyear=dyear, year_range=year_range,
         years_ch=years_ch, years=years)

}
