#' Is species coastal
#'
#' @description This function is used to determine whether a species is coastal.
#' It fits a simple occupancy model with distance to coast as a covariate for
#' occupancy probability. If the coefficient is large (>2) then the species is
#' considered coastal and a diagnostic file is created in a directory selected
#' by the user.
#' @param sp_sel The SAFRING code of the species being assessed.
#' @param out_dir The output directory where information must be saved and where
#' the diagnostic file should be found.
#' @param visit_mod The model for the detection probability. The occupancy model
#' is fixed to `psi ~ 1 + dist_coast`.
#' @param force Logical; if FALSE the function will check if the an assessment
#' has been done for this species already and skip if it has. See details.
#'
#' @return It returns TRUE if the species is coastal and it will create a file
#' in `out_dir` to indicate that the check has already been done for this
#' species. In the future, the function will check if this file exits before
#' and avoid fitting the model if it does, unless force == TRUE
#' @export
#'
#' @examples
#' \dontrun{
#' isSpCoastal(235, outdir, reformulate(c("1", "log(TotalHours+1)"), response = "p"))
#' }
isSpCoastal <- function(sp_sel, out_dir, visit_mod){

    coast_idx <- grep("^coast", list.files(paste0(out_dir, sp_sel)))
    coast_file <- list.files(paste0(out_dir, sp_sel))[coast_idx]

    if(length(coast_idx) == 0){

        # Fit a simple model with distance to coast as a covariate
        print("Checking if species is coastal...")

        m1 <- fit_occu(forms = list(visit_mod,
                                    reformulate(c("1", "dist_coast"), response = "psi")),
                       visit_data = occuRdata$visit,
                       site_data = occuRdata$site,
                       print = FALSE)

        # If coefficient for dist_coast is really large then spp is coastal
        if(m1$fit$par[2] < -2){
            out <- TRUE
            file.create(paste0(out_dir, sp_sel, "/coast_", sp_sel,"_", "1.txt"))
        } else {
            out <- FALSE
            file.create(paste0(out_dir, sp_sel, "/coast_", sp_sel,"_", "0.txt"))
        }

    } else if(coast_file == paste0("coast_", sp_sel,"_", "1.txt")){
        out <- TRUE
    } else if(coast_file == paste0("coast_", sp_sel,"_", "0.txt")){
        out <- FALSE
    }

    return(out)

}
