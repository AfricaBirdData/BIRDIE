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
#' @param ncores Number of cores to be used for the extraction. Defaults to 1
#'
#' @return The sites dataframe with an additional column corresponding to
#' the new variable.
#' @export
#'
#' @examples
addOccVisitCovt <- function(visits, sites, covt, covts_dir, file_fix,
                            ncores = 1){

    # Identify all different year-month combinations
    visits <- visits %>%
        dplyr::mutate(month = dplyr::if_else(month < 10,
                                             paste0("0", month),
                                             as.character(month))) %>%
        dplyr::mutate(yearmonth = paste0("X", year, ".", month))

    yearmonth <- visits %>%
        dplyr::distinct(yearmonth) %>%
        dplyr::pull(yearmonth)

    covt_r <- readRDS(paste0(covts_dir, file_fix[1], covt, file_fix[2], ".rds"))

    covt_r <- raster::subset(covt_r,
                             grep(paste(yearmonth, collapse = "|"),
                                  names(covt_r), value = TRUE))

    # Need to subset sites for each year month
    visits <- dplyr::nest_by(visits, yearmonth)

    # Transform to spatial sites to sp package to avoid problems with furrr
    sites <- sf::as_Spatial(sites)

    print(paste("Extracting variable", covt))

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

    #for(i in seq_len(nrow(visits))){

    f <- function(.yr, .visits, .sites=sites, .covt_r=covt_r, .covt=covt){
        # remove the next line
        # .yr=visits$yearmonth[4];
        # .visits=visits$data[[4]]
        # .sites=sites
        # .covt_r=covt_r
        # .covt=covt

        ry <- raster::subset(.covt_r, grep(.yr, names(.covt_r)))

        subsites <- .sites[.sites$Name %in% .visits$SiteName,]

        vv <- raster::extract(ry, subsites, fun = mean)

        subsites <- as.data.frame(subsites) %>%
            dplyr::mutate(!!.covt := vv) %>%
            dplyr::select(SiteName = Name, dplyr::all_of(.covt))

        covt_out <- .visits %>%
            dplyr::left_join(subsites, by = "SiteName")

        return(covt_out)
    }

    out <- furrr::future_map2_dfr(visits$yearmonth,
                                  visits$data,
                                  ~f(.yr = .x, .visits = .y))

    return(out)
}
