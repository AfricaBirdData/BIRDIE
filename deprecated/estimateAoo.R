#' Estimate Area of Occupancy (AOO)
#'
#' @param year Year selected to estimate AOO for.
#' @param spp_code SAFING_No for the species predicting AOO for, as coming from
#' the ABAP API.
#' @param site_data Site data used to fit the occupancy model from which AOO is
#' estimated.
#' @param visit_data Visit data used to fit the occupancy model from which AOO is
#' estimated
#' @param config A named list coming from \link{configPreambOccuR}.
#' @param verbose To silence data update message, set to FALSE.
#'
#' @return The indicator file will be updated with the resulting AOO info, which
#' includes: spp_code, year, estimate, sd, lower 95% CI, upper 95% CI
#' @export
#'
#' @examples
estimateAoo <- function(year, spp_code, site_data, visit_data, config,
                        verbose = TRUE){

    # Format data to occuR
    occuRdata <- BIRDIE::prepDataOccuR(spp_code = spp_code,
                                       years = config$years,
                                       site_data = site_data,
                                       visit_data = visit_data,
                                       download = TRUE)

    # Load fit
    fit <- with(config,
                readRDS(file.path(fit_dir, spp_code, paste0("occur_fit_", years_ch, "_", spp_code, ".rds"))))

    # Prepare prediction data
    pred_data <- BIRDIE::prepPredictDataOccuR(occuRdata,
                                              sf::st_drop_geometry(site_data),
                                              years = year,
                                              scaling = TRUE)

    # Predict
    pred <- predict(fit, occuRdata$visit, pred_data, nboot = 1000)


    # Estimate realized occupancy ---------------------------------------------

    occu <- matrix(nrow = 1000, ncol = nrow(site_data))

    for(i in 1:1000){

        pred_occu <- BIRDIE::summarizePredOccuR(pred_p = pred$pboot[i,],
                                                pred_psi = pred$psiboot[i,],
                                                pred_data = pred_data,
                                                visit_data = occuRdata$visit)

        occu[i,] <- pred_occu$real_occu

    }

    aoo <- rowSums(occu)

    # retrieve indicator file
    ind_file <- read.csv(file.path(config$fit_dir, spp_code, paste0("indtr_", spp_code,".csv")))

    # Add new row
    ind_file[nrow(ind_file) + 1, ] <- c(spp_code,
                                        "aoo",
                                        paste("01", "01", year, sep = "-"),
                                        paste("31", "12", year, sep = "-"),
                                        round(mean(aoo), 3),
                                        round(sd(aoo), 3),
                                        round(quantile(aoo, 0.025), 3),
                                        round(quantile(aoo, 0.975), 3))

    write.csv(ind_file,
              file.path(config$fit_dir, spp_code, paste0("indtr_", spp_code,".csv")),
              row.names = FALSE)

    # Print results
    if(verbose){
        print(paste(file.path(config$fit_dir, spp_code, paste0("indtr_", spp_code,".csv")), "updated with"))
        print(ind_file[nrow(ind_file), ])
    }


}
