#' Fast raster extraction in parallel
#'
#' @description This function is a wrapper around exactextractr::exactextract()
#' to perform extraction of raster values within polygons in parallel. It
#' divides the data into blocks and performs extraction using one core per block.
#' @param rst A raster or raster stack/brick
#' @param spt A spatial polygon object (sf, sfc, SpatialPolygonDataFrame)
#' @param ncores Number of cores to use in the extraction
#' @param fun Character string naming the function to use for summarizing the
#' raster values per polygon.
#' @param ... Other arguments to pass onto the function 'fun'.
#' exactextractr::exactextract()
#'
#' @return An spatial polygon object with an extract column containing the
#' summary of raster values per polygon
#' @export
#'
#' @examples
exactExtractParll <- function(rst, spt, ncores, fun, ... ){

    if(!"raster" %in% .packages()){
        stop("Package raster needs to be loaded")
    }

    spt$grp <- ceiling(seq_len(nrow(spt)) / (nrow(spt) / ncores))

    spt <- spt %>%
        dplyr::group_split(grp)

    vals <- furrr::future_map_dfr(spt,
                                  ~exactextractr::exact_extract(x = rst,
                                                                y = .x,
                                                                fun = fun,
                                                                append_cols = TRUE),
                                  .options = furrr::furrr_options(seed = TRUE,
                                                                  packages = c("raster")))

    out <- spt %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(vals = vals[,fun]) %>%
        dplyr::select(-grp)

    return(out)
}
