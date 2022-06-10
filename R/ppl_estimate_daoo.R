#' Estimate change in Area of Occupancy
#'
#' @inheritParams ppl_run_pipe_dst2
#' @param term A character string with one of three options with the time-frame
#' for the change: 'annual' (1 year), 'short' (15 years), 'long' (30 years).
#' @param verbose Logical. If TRUE, the new line added to the indicator table
#' is displayed in the console.
#' @param overwrite Logical, If TRUE, previous DAOO calculations will be overwritten.
#'
#' @return
#' @export
#'
#' @examples
ppl_estimate_daoo <- function(sp_code, config, term = c("annual", "short", "long"),
                              verbose, overwrite, ...){

    dyear <- dplyr::case_when(term == "annual" ~ 1,
                              term == "short" ~ 10, # This should be 15 but still have no data in SABAP2 for 15 years
                              term == "long" ~ 30)

    # retrieve indicator files
    indtr_file1 <- file.path(config$out_dir, sp_code, paste0("indtr_dst_", sp_code, "_", config$year - dyear, ".csv"))
    indtr_file2 <- file.path(config$out_dir, sp_code, paste0("indtr_dst_", sp_code, "_", config$year, ".csv"))

    if(file.exists(indtr_file1)){
        indtr <- rbind(utils::read.csv(indtr_file1),
                       utils::read.csv(indtr_file2))
    } else {
        indtr <- utils::read.csv(indtr_file2)
    }

    if(overwrite){
        indtr <- indtr %>%
            dplyr::filter(indicator != "daoo")
    }

    new_rows <- vector("list", length(config$years))

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
        daoo <- daoo %>%
            dplyr::group_by(start_date, end_date) %>%
            dplyr::filter(opt == min(opt)) %>%
            dplyr::ungroup()

        opt_daoo <- max(daoo$opt)

        daoo <- daoo %>%
            dplyr::mutate(start_date = as.Date(start_date, format = "%d-%m-%Y")) %>%
            dplyr::arrange(start_date) %>%
            dplyr::mutate(estimate = estimate - dplyr::lag(estimate),
                             st_dev = sqrt(st_dev^2 + dplyr::lag(st_dev^2))) %>%
            dplyr::filter(!is.na(estimate)) %>%
            dplyr::mutate(indicator = "daoo",
                          term = term,
                          start_date = ini,
                          end_date = end,
                          estimate = round(estimate, 2),
                          st_dev = round(st_dev, 3),
                          lb95 = round(qnorm(0.025, estimate, st_dev), 3),
                          ub95 = round(qnorm(0.975, estimate, st_dev), 3),
                          opt = opt_daoo)

        new_rows[[t]] <- daoo

        # Print results
        if(verbose){
            message(paste(indtr_file2, "updated with"))
            message(new_rows[[t]])
        }
    }

    # Add new rows to last indicator file
    new_rows <- dplyr::bind_rows(new_rows)

    out <- utils::read.csv(indtr_file2)

    if(overwrite){
        out <- out %>%
            dplyr::filter(indicator != "daoo")
    }

    rbind(out, new_rows) %>%
        dplyr::arrange(indicator, term, start_date) %>%
        utils::write.csv(indtr_file2,
                         row.names = FALSE)

}
