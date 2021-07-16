#' Add occupancy site covariate
#'
#' @param sites A spatial object with the sites where data was collected.
#' @param covt A character string with the name of the covariate to be
#' extracted. The names must correspond to the names contained the covariate
#' file name (e.g. "prcp", "tmax", etc.). At the moment the function assumes
#' that covariates are found at: "covts_dir_file_fix[1]covt_file_fix[2].rds"
#' @param covts_dir Directory where covariate are found.
#' @param years A vector of years that we want to extract
#' @param file_fix A character vector with two elements corresponding to
#' additional characters that are found before and after covts.
#' @param ncores Number of cores to be used for the extraction. Defaults to 1
#'
#' @return The sites dataframe with an additional column corresponding to
#' the new variable.
#' @export
#'
#' @examples
addOccSiteCovt <- function(sites, covt, years, covts_dir, file_fix, ncores = 1){

    covt_r <- readRDS(paste0(covts_dir, file_fix[1], covt,
                           file_fix[2], ".rds"))

    covt_r <- raster::subset(covt_r,
                             grep(paste(years, collapse = "|"),
                                  names(covt_r), value = TRUE))

    # Set up multicore processing
    if (ncores > 1) {
        pnow <- future::plan()
        on.exit(future::plan(pnow), add = TRUE)

        if (ncores > future::availableCores()) {
            stop("Not enough available cores")
        } else {
            future::plan("multiprocess", workers = ncores)
        }
    }

    print(paste("Extracting variable", covt))

    f <- function(.yr, .covt_r=covt_r, .covt=covt, .sites=sites){

        ry <- raster::subset(.covt_r, grep(.yr, names(.covt_r)))
        ry <- raster::calc(ry, fun = mean)
        names(ry) <- .yr

        covt_name <- paste0(.covt, "_", .yr)

        covt_out <- .sites %>%
            sf::st_drop_geometry() %>%
            dplyr::mutate(!!covt_name := raster::extract(ry, as(.sites, "Spatial"), fun = mean)) %>%
            dplyr::select(dplyr::all_of(covt_name))

        return(covt_out)
    }

    out <- furrr::future_map_dfc(years, ~f(.x))

    out <- cbind(sites, out)

    return(out)
}
