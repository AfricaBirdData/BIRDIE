#' Estimate change in Area of Occupancy
#'
#' @inheritParams ppl_run_pipe_dst2
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
    indtr_file <- file.path(config$out_dir, sp_code, paste0("indtr_dst_", sp_code, "_", year, ".csv"))
    indtr <- utils::read.csv(indtr_file)

    for(t in seq_along(config$years)){

        year <- config$years[t]

        if(year < 2009){
            next
        }

        # Set initial date
        ini_year <- dplyr::case_when(term == "annual" ~ year - 1,
                                     term == "short" ~ year - 10, # This should be 15 but still have no data in SABAP2 for 15 years
                                     term == "long" ~ year - 30)

        ini <- paste("01", "01", ini_year, sep = "-")
        end <- paste("31", "12", year, sep = "-")

        # Estimate daoo
        daoo <- indtr %>%
            dplyr::filter(species == sp_code & indicator == "aoo", (start_date == ini | end_date == end))

        if(nrow(daoo) < 2){
            print(paste("We are probably missing AOO for previous years to compute dAOO in", year))
            next
        }

        # Set opt value to max value found for AOO
        opt_daoo <- max(daoo$opt)

        daoo <- daoo %>%
            dplyr::mutate(start_date = as.Date(start_date, format = "%d-%m-%Y")) %>%
            dplyr::arrange(start_date) %>%
            dplyr::mutate(estimate = estimate - dplyr::lag(estimate),
                             st_dev = sqrt(st_dev^2 + dplyr::lag(st_dev^2))) %>%
            dplyr::filter(!is.na(estimate)) %>%
            dplyr::mutate(indicator = "daoo",
                          term = tt,
                          start_date = ini,
                          end_date = end,
                          estimate = round(estimate, 2),
                          st_dev = round(st_dev, 3),
                          lb95 = round(qnorm(0.025, estimate, st_dev), 3),
                          ub95 = round(qnorm(0.975, estimate, st_dev), 3),
                          opt = opt_daoo)

        # Add new row
        indtr <- rbind(indtr, daoo)

        # Print results
        if(verbose){
            print(paste(indtr_file, "updated with"))
            print(indtr[nrow(indtr), ])
        }
    }

    indtr %>%
        dplyr::arrange(indicator, term, start_date) %>%
        utils::write.csv(indtr_file,
                         row.names = FALSE)

}
