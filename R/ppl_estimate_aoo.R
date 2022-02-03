#' Estimate Area of Occupancy
#'
#' @inheritParams ppl_run_pipe_distr
#' @param verbose Logical. If TRUE, the new line added to the indicator table
#' is displayed in the console.
#'
#' @return
#' @export
#'
#' @examples
ppl_estimate_aoo <- function(sp_code, year, config, verbose, ...){

    # Calculate the position of the estimated year in relation to the optimal year
    # This will be used for the opt field
    opt <- abs(year - config$year)

    # Set initial date and final date
    ini <- paste("01", "01", year, sep = "-")
    end <- paste("31", "12", year, sep = "-")

    # retrieve indicator file
    ind_file <- read.csv(file.path(config$fit_dir, sp_code, paste0("indtr_", sp_code,".csv")))

    # Check if the record already exists
    case <- ind_file %>%
        dplyr::filter(species == sp_code & indicator == "aoo" &
                          start_date == ini & end_date == end)

    # If it exists and current opt is smaller than the new opt stop
    if(nrow(case) > 0 && case$opt <= opt){
        return(warning(paste("Year", year, "AOO not updated because there is a better record in the database")))
    } else if(nrow(case) > 0 && case$opt > opt){
        ind_file <- ind_file %>%
            dplyr::filter(!(species == sp_code & indicator == "aoo" &
                                start_date == ini & end_date == end))
    }

    # Predict from model
    preds <- ppl_predict_occur(sp_code, year, config, ...)


    # Estimate realized occupancy ---------------------------------------------

    print("Estimating Area of Occupancy")

    occu <- matrix(nrow = 1000, ncol = nrow(preds$pred_data))

    for(i in 1:1000){

        pred_occu <- summarizePredOccuR(pred_p = preds$pred$pboot[i,],
                                        pred_psi = preds$pred$psiboot[i,],
                                        pred_data = preds$pred_data,
                                        visit_data = preds$occuRdata$visit)

        occu[i,] <- pred_occu$real_occu

    }

    aoo <- rowSums(occu)

    # Add new row
    ind_file[nrow(ind_file) + 1, ] <- c(sp_code,
                                        "aoo",
                                        ini,
                                        end,
                                        "annual",
                                        round(mean(aoo), 3),
                                        round(sd(aoo), 3),
                                        round(quantile(aoo, 0.025), 3),
                                        round(quantile(aoo, 0.975), 3),
                                        opt)

    write.csv(ind_file,
              file.path(config$fit_dir, sp_code, paste0("indtr_", sp_code,".csv")),
              row.names = FALSE)

    # Print results
    if(verbose){
        print(paste(file.path(config$fit_dir, sp_code, paste0("indtr_", sp_code,".csv")), "updated with"))
        print(ind_file[nrow(ind_file), ])
    }

}
