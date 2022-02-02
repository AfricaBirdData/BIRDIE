#' Estimate change in Area of Occupancy
#'
#' @inheritParams ppl_run_pipe_distr
#' @param term A character string with one of three options with the time-frame
#' for the change: 'annual' (1 year), 'short' (15 years), 'long' (30 years).
#' @param verbose Logical. If TRUE, the new line added to the indicator table
#' is displayed in the console.
#'
#' @return
#' @export
#'
#' @examples
ppl_estimate_daoo <- function(sp_code, year, config, term = c("annual", "short", "long"), verbose, ...){

    # Set initial date
    ini_year <- dplyr::case_when(term == "annual" ~ year - 1,
                                 term == "short" ~ year - 10, # This should be 15 but still have no data in SABAP2 for 15 years
                                 term == "long" ~ year - 30)

    ini <- paste("01", "01", ini_year, sep = "-")
    end <- paste("31", "12", year, sep = "-")

    # Calculate the position of the estimated year in relation to the optimal year
    # This will be used for the opt field
    opt <- abs(year - config$year)

    # retrieve indicator file
    ind_file <- read.csv(file.path(config$fit_dir, sp_code, paste0("indtr_", sp_code,".csv")))

    # Check if the record already exists
    tt <- term # This is just to ensure stability
    case <- ind_file %>%
        dplyr::filter(species == sp_code & indicator == "daoo" &
                          start_date == ini & end_date == end & term == tt)

    # If it exists and current opt is smaller than the new opt stop
    if(nrow(case) > 0 && case$opt <= opt){
        return(warning(paste("Year", year, "DAOO not updated because there is a better record in the database")))
    } else if(nrow(case) > 0 && case$opt > opt){
        ind_file <- ind_file %>%
            dplyr::filter(!(species == sp_code & indicator == "daoo" &
                                start_date == ini & end_date == end & term == tt))
    }

    # Estimate daoo
    daoo <- ind_file %>%
        dplyr::filter(species == sp_code & indicator == "aoo", (start_date == ini | end_date == end)) %>%
        dplyr::mutate(start_date = as.Date(start_date, format = "%d-%m-%Y")) %>%
        dplyr::arrange(start_date) %>%
        dplyr::transmute(estimate = estimate - dplyr::lag(estimate),
                         st_dev = sqrt(st_dev^2 + dplyr::lag(st_dev^2))) %>%
        dplyr::filter(!is.na(estimate))

    # Add new row
    ind_file[nrow(ind_file) + 1, ] <- c(sp_code,
                                        "daoo",
                                        ini,
                                        end,
                                        term = term,
                                        round(daoo$estimate, 3),
                                        round(daoo$st_dev, 3),
                                        round(qnorm(0.025, daoo$estimate, daoo$st_dev), 3),
                                        round(qnorm(0.975, daoo$estimate, daoo$st_dev), 3),
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
