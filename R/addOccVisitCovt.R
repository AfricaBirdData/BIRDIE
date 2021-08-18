#' Add occupancy visit covariate
#'
#' @param visits A dataframe with detection/non-detection data. It must contain
#' at least the following columns: "year", "month", "SiteName".
#' @param sites A spatial object with the sites where data was collected.
#' @param covt A character string with the names of the covariate to be
#' extracted. The names must correspond to the names contained the covariate
#' file name (e.g. "prcp", "tmax", etc.). At the moment the function assumes
#' that covariates are found at: "covts_dir_file_fix[1]covt_file_fix[2].rds"
#' @param covts_dir Directory where covariate are found.
#' @param file_fix A character vector with two elements corresponding to
#' additional characters that are found before and after covts.
#'
#' @return The sites dataframe with an additional column corresponding to
#' the new variable.
#' @export
#'
#' @examples
addOccVisitCovt <- function(visits, sites, covt, covts_dir, file_fix){

    # Identify all different year-month combinations
    visits <- visits %>%
        dplyr::mutate(month_p = dplyr::if_else(month < 10,
                                             paste0("0", month),
                                             as.character(month))) %>%
        dplyr::mutate(yearmonth = paste0("X", year, ".", month_p))

    yearmonth <- visits %>%
        dplyr::distinct(yearmonth) %>%
        dplyr::pull(yearmonth)

    covt_r <- readRDS(paste0(covts_dir, file_fix[1], covt, file_fix[2], ".rds"))

    covt_r <- raster::subset(covt_r,
                             grep(paste(yearmonth, collapse = "|"),
                                  names(covt_r), value = TRUE))

    # Need to subset sites for each year month
    visits <- dplyr::nest_by(visits, yearmonth)

    print(paste("Extracting variable", covt))

    f <- function(.yr, .visits, .sites=sites, .covt_r=covt_r, .covt=covt){

        ry <- raster::subset(.covt_r, grep(.yr, names(.covt_r)))

        subsites <- .sites[.sites$Name %in% .visits$SiteName,]

        vv <- exactextractr::exact_extract(ry, subsites, fun = "mean",
                                           progress = FALSE)

        subsites <- as.data.frame(subsites) %>%
            dplyr::mutate(!!.covt := vv) %>%
            dplyr::select(SiteName = Name, dplyr::all_of(.covt))

        covt_out <- .visits %>%
            dplyr::left_join(subsites, by = "SiteName")

        return(covt_out)
    }

    out <- furrr::future_map2_dfr(visits$yearmonth,
                                  visits$data,
                                  ~f(.yr = .x, .visits = .y),
                                  .options = furrr::furrr_options(seed = TRUE,
                                                                  packages = c("sf", "raster")))

    return(out)
}
