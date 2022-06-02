#' Estimate Area of Occupancy
#'
#' @inheritParams ppl_run_pipe_dst2
#' @param verbose Logical. If TRUE, the new line added to the indicator table
#' is displayed in the console.
#'
#' @return
#' @export
#'
#' @examples
ppl_estimate_aoo <- function(sp_code, config, verbose, ...){

    # retrieve indicator file
    indtr_file <- file.path(config$out_dir, sp_code, paste0("indtr_dst_", sp_code, "_", year, ".csv"))
    indtr <- utils::read.csv(indtr_file)

    # Calculate the position of year in relation to the optimal year
    # This will be used for the opt field
    opts <- abs(1:config$dur - ceiling(config$dur/2))

    for(t in seq_along(config$years)){

        year <- config$years[t]

        opt <- opts[t]

        # Set initial date and final date
        ini <- paste("01", "01", year, sep = "-")
        end <- paste("31", "12", year, sep = "-")

        # Check if the record already exists
        case <- indtr %>%
            dplyr::filter(species == sp_code & indicator == "aoo" &
                              start_date == ini & end_date == end)

        # If it exists and current opt is smaller than the new opt stop
        if(nrow(case) > 0 && case$opt <= opt){
            print(paste("Year", year, "AOO not updated because there is a better record in the database"))
            next
        } else if(nrow(case) > 0 && case$opt > opt){
            indtr <- indtr %>%
                dplyr::filter(!(species == sp_code & indicator == "aoo" &
                                    start_date == ini & end_date == end))
        }

        # Predict from model
        preds <- predictOccuR(sp_code, year, config, ...)

        # Stop if there are too few detections
        if(is.numeric(preds) && preds %in% c(1, 2)){
            return(1)
        }


        # Estimate realized occupancy ---------------------------------------------

        print("Estimating Area of Occupancy")

        occu <- matrix(nrow = 1000, ncol = nrow(preds$pred_data))

        for(i in 1:1000){

            pred_occu <- summarizeOccuRPred(pred_p = preds$pred$pboot[i,],
                                            pred_psi = preds$pred$psiboot[i,],
                                            pred_data = preds$pred_data,
                                            visit_data = preds$occuRdata$visit)

            occu[i,] <- pred_occu$real_occu

        }

        aoo <- rowSums(occu)

        # Add new row
        new_row <- data.frame(species = sp_code,
                              indicator = "aoo",
                              start_date = ini,
                              end_date = end,
                              term = "annual",
                              estimate = round(mean(aoo), 3),
                              st_dev = round(sd(aoo), 3),
                              lb95 = round(quantile(aoo, 0.025), 3),
                              ub95 = round(quantile(aoo, 0.975), 3),
                              opt = opt)

        indtr <- rbind(indtr, new_row)

        # Print results
        if(verbose){
            print(paste(indtr_file, "updated with"))
            print(indtr[nrow(indtr), ])
        }

    }

    # Save to disc
    indtr %>%
        dplyr::arrange(indicator, term, start_date) %>%
        utils::write.csv(indtr_file,
                         row.names = FALSE)

    return(0)

}
