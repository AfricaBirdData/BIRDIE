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
ppl_estimate_daoo <- function(sp_code, config, term = c("annual", "short", "long"), verbose, ...){

    # retrieve indicator file
    indtr_file <- file.path(config$out_dir, sp_code, paste0("indtr_dst_", sp_code,".csv"))
    indtr <- read.csv(indtr_file)

    # Calculate the position of year in relation to the optimal year
    # This will be used for the opt field
    opts <- abs(1:config$dur - ceiling(config$dur/2))

    for(t in seq_along(config$years)){

        year <- config$years[t]

        if(year < 2009){
            next
        }

        opt <- opts[t]

        # Set initial date
        ini_year <- dplyr::case_when(term == "annual" ~ year - 1,
                                     term == "short" ~ year - 10, # This should be 15 but still have no data in SABAP2 for 15 years
                                     term == "long" ~ year - 30)

        ini <- paste("01", "01", ini_year, sep = "-")
        end <- paste("31", "12", year, sep = "-")

        # Check if the record already exists
        tt <- term # This is just to ensure stability
        case <- indtr %>%
            dplyr::filter(species == sp_code & indicator == "daoo" &
                              start_date == ini & end_date == end & term == tt)

        # If it exists and current opt is smaller than the new opt stop
        if(nrow(case) > 0 && case$opt <= opt){
            return(warning(paste("Year", year, "DAOO not updated because there is a better record in the database")))
        } else if(nrow(case) > 0 && case$opt > opt){
            indtr <- indtr %>%
                dplyr::filter(!(species == sp_code & indicator == "daoo" &
                                    start_date == ini & end_date == end & term == tt))
        }

        # Estimate daoo
        daoo <- indtr %>%
            dplyr::filter(species == sp_code & indicator == "aoo", (start_date == ini | end_date == end)) %>%
            dplyr::mutate(start_date = as.Date(start_date, format = "%d-%m-%Y")) %>%
            dplyr::arrange(start_date) %>%
            dplyr::mutate(estimate = estimate - dplyr::lag(estimate),
                             st_dev = sqrt(st_dev^2 + dplyr::lag(st_dev^2))) %>%
            dplyr::filter(!is.na(estimate)) %>%
            dplyr::mutate(indicator = "daoo",
                          term = tt,
                          start_date = as.character(start_date, format = "%d-%m-%Y"),
                          estimate = round(estimate, 2),
                          st_dev = round(st_dev, 3),
                          lb95 = round(qnorm(0.025, estimate, st_dev), 3),
                          ub95 = round(qnorm(0.975, estimate, st_dev), 3),
                          opt = opt)

        if(nrow(daoo) == 0){
            warning(paste("We are probably missing AOO for previous years to compute dAOO in", year))
            next
        }

        # Add new row
        indtr <- rbind(indtr, daoo)

        # Print results
        if(verbose){
            print(paste(indtr_file, "updated with"))
            print(indtr[nrow(indtr), ])
        }
    }

    write.csv(indtr,
              indtr_file,
              row.names = FALSE)

}
