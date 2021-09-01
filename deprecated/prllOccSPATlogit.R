#' Parallel Spatial Occupancy with Rcppocc
#'
#' @description Wrapper around Rcppocc::occSPATlogit to run multiple chains.
#' They may run in parallel if future::plan is configured to do so previously.
#' @param nchains Number of chains to run.
#' @param ... Arguments to be used by Rcppocc::occSPATlogit
#'
#' @return A list with nchains outputs from Rcppocc
#' @export
#'
#' @examples
prllOccSPATlogit <- function(nchains, ...){

    samples <- furrr::future_map(seq_len(nchains), function(.x,...) occSPATlogit(...), ...,
                             .options = furrr::furrr_options(seed = TRUE))

    npar <- length(samples[[1]])

    out <- vector("list", length = npar)

    for(i in seq_len(npar)){

        a <- lapply(samples, "[[", i)

        out[[i]] <- array(dim = c(nrow(a[[1]]), ncol(a[[1]]), length(a)))

        for(j in seq_along(a)){
            out[[i]][,,j] <- a[[j]]
        }

    }

    names(out) <- names(samples[[1]])

    return(out)

}
