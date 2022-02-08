#' Estimate Area Of Occupancy for species type
#'
#' @inheritParams ppl_estimate_distr_sp_type
#'
#' @return
#'
#' @examples
estimateAOOSpType <- function(sp_type, year, config, force_predict = FALSE, verbose){

    # Lookup species in the group
    if(sp_type == "group"){
        species <- c(4, 6, 41, 235)
    } else {
        stop("Sorry the group doesn't exist")
    }

    # Calculate the position of the estimated year in relation to the optimal year
    # This will be used for the opt field
    opt <- abs(year - config$year)

    # Set initial date and final date
    ini <- paste("01", "01", year, sep = "-")
    end <- paste("31", "12", year, sep = "-")

    # retrieve indicator file
    ind_file <- read.csv(file.path(config$fit_dir, sp_type, paste0("indtr_", sp_type,".csv")))

    # Check if the record already exists
    case <- ind_file %>%
        dplyr::filter(species == sp_type & indicator == "aoo" &
                          start_date == ini & end_date == end)

    # If it exists and current opt is smaller than the new opt stop
    if(nrow(case) > 0 && case$opt <= opt){
        return(warning(paste("Year", year, "AOO not updated because there is a better record in the database")))
    } else if(nrow(case) > 0 && case$opt > opt){
        ind_file <- ind_file %>%
            dplyr::filter(!(species == sp_type & indicator == "aoo" &
                                start_date == ini & end_date == end))
    }

    # Loop through group species
    for(s in seq_along(species)){

        sp_code <- species[s]

        if(!file.exists(file.path(tempdir(), paste0(sp_code, "_", year, "_occu.rds"))) | force_predict){

            # Predict from model
            preds <- ppl_predict_occur(sp_code, year, config)

            # Estimate realized occupancy

            print(paste("Estimating Area of Occupancy for species", sp_code))

            occu <- matrix(nrow = 1000, ncol = nrow(preds$pred_data))

            for(i in 1:1000){

                pred_occu <- summarizePredOccuR(pred_p = preds$pred$pboot[i,],
                                                pred_psi = preds$pred$psiboot[i,],
                                                pred_data = preds$pred_data,
                                                visit_data = preds$occuRdata$visit)

                occu[i,] <- pred_occu$real_occu

            }

            saveRDS(occu, file.path(tempdir(), paste0(sp_code, "_", year, "_occu.rds")))
        }

    }

    # Find probability of not detecting any of the species and find reciprocal (detect at least one)
    preds <- vector("list", length = length(species))

    for(s in seq_along(species)){

        preds[[s]] <- 1 - readRDS(file.path(tempdir(), paste0(species[s], "_", year, "_occu.rds")))

    }

    preds <- Reduce("*", preds)

    occu <- 1 - preds

    aoo <- rowSums(occu)

    # Add new row
    ind_file[nrow(ind_file) + 1, ] <- c(sp_type,
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
              file.path(config$fit_dir, sp_type, paste0("indtr_", sp_type,".csv")),
              row.names = FALSE)

    # Print results
    if(verbose){
        print(paste(file.path(config$fit_dir, sp_type, paste0("indtr_", sp_type,".csv")), "updated with"))
        print(ind_file[nrow(ind_file), ])
    }

}
