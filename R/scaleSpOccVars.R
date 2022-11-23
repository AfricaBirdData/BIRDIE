#' Scale covariates in spOccupancy-type data
#'
#' @param spOcc_data an spOccupancy data list.
#' @param var_type Type of variables we want to scale. Currently, one of "occ"
#' occupancy covariates, "det" detection covariates.
#' @param scale_vars A vector with the names of the covariates that we want to
#' scale.
#'
#' @return An spOccupancy data list with the scaled covariates, substituting the
#' original, unscaled covariates. A single factor is used to center and scale all
#' data (across all dimensions) of the covariate. That means, for example, that
#' all seasons will be scaled by the same amount. The factors used
#' for centering and scaling are stored as attributes of each covariate.
#' @export
#'
#' @examples
scaleSpOccVars <- function(spOcc_data, var_type, scale_vars){

    # if(var_type == "det"){
    #     stop("Only var_type = 'occ' is supported at the moment")
    # }

    if(var_type == "occ"){

        f <- function(x){
            scl <- stats::sd(c(x), na.rm = TRUE)
            cnt <- mean(c(x), na.rm = TRUE)

            out_scl <- apply(x, 2, FUN = scale, center = cnt, scale = scl)
            dimnames(out_scl) <- dimnames(x)
            attr(out_scl, "scaled:scale") <- scl
            attr(out_scl, "scaled:center") <- cnt
            out_scl
        }

        covts <- spOcc_data$occ.covs[scale_vars]

        covts <- lapply(covts, f)

        spOcc_data$occ.covs[scale_vars] <- covts

        } else if(var_type == "det"){

            f <- function(x){
                scl <- stats::sd(c(x), na.rm = TRUE)
                cnt <- mean(c(x), na.rm = TRUE)

                out_scl <- apply(x, 2, FUN = scale, center = cnt, scale = scl)
                dimnames(out_scl) <- dimnames(x)
                attr(out_scl, "scaled:scale") <- scl
                attr(out_scl, "scaled:center") <- cnt
                out_scl
            }

            covts <- spOcc_data$det.covs[scale_vars]

            covts <- lapply(covts, f)

            spOcc_data$det.covs[scale_vars] <- covts

        }

    return(spOcc_data)

    }
