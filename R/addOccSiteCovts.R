#' Add occupancy site covariate
#'
#' @param sites A spatial object with the sites where data was collected.
#' @param covt A character string with the names of the covariate to be
#' extracted. The names must correspond to the names contained the covariate
#' file name (e.g. "prcp", "tmax", etc.). At the moment the function assumes
#' that covariates are found at: "covts_dir_file_fix[1]covt_file_fix[2].rds"
#' @param covts_dir Directory where covariate are found.
#' @param years A vector of years that we want to extract
#' @param file_fix A character vector with two elements corresponding to
#' additional characters that are found before and after covts.
#'
#' @returnThe sites dataframe with additional columns corresponding to
#' the new variables.
#' @export
#'
#' @examples
addOccSiteCovt <- function(sites, covt, years, covts_dir, file_fix){

    covt_r <- readRDS(paste0(covts_dir, file_fix[1], covt,
                           file_fix[2], ".rds"))

    covt_r <- raster::subset(covt_r,
                             grep(paste(years, collapse = "|"),
                                  names(covt_r), value = TRUE))

    # covt_r <- raster::subset(.covt_r, grep(years, names(.covt_r)))
    # yrs <- raster::nlayers(covt_r)/12

    print(paste("Extracting variable", covt))

    f <- function(.yr, .covt_r=covt_r, .covt=covt, .sites=sites){

        ry <- raster::subset(.covt_r, grep(.yr, names(.covt_r)))
        ry <- raster::calc(ry, fun = mean)
        names(ry) <- .yr
        rr <- ry

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
