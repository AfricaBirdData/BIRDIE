#' Prepare covariate data for abundance pipeline
#'
#' @inheritParams ppl_run_pipe_abu1
#' @param catchment A sf object with polygons corresponding to the catchment, or
#' reference area considered for the covariates associated with CWAC sites.
#' @param ... Other arguments passed on to \link{prepGEESpCountData}
#'
#' @return
#' @export
#'
#' @examples
ppl_create_data_ssm <- function(sp_code, year, catchment, config, ...){

    # Download species data
    print("Downloading from CWAC")

    sp_data <- CWAC::getCwacSppCounts(sp_code)

    # In what South African sites was the species present?
    sp_sites <- sp_data %>%
        dplyr::filter(Country == "South Africa") %>%
        dplyr::pull(LocationCode) %>%
        unique()


    # Add missing counts ------------------------------------------------------

    print("Adding missing counts")

    # List to store data from multiple sites
    counts <- vector("list", length = length(sp_sites))

    # For each site:
    for(s in seq_along(sp_sites)){

        site_sel <- sp_sites[s]

        # Temp file to save to/retrieve from
        temp_file <- file.path(tempdir(), paste(site_sel, config$years_ch, "cwac_dat.rds", sep = "_"))

        if(file.exists(temp_file)){

            site_data_w_miss <- readRDS(temp_file)

        } else {

            # Download site data and prepare missing visits
            site_data <- CWAC::getCwacSiteCounts(site_sel)

            # Figure out which counts are missing and which are zero for the species of interest
            site_data_w_miss <- CWAC::addMissingCwacCounts(site_data, config$years)

            # Give missing surveys a date based on the dates from other surveys
            month_summer <- site_data %>%
                dplyr::mutate(month = lubridate::month(StartDate)) %>%
                dplyr::filter(Season == "S", !is.na(month)) %>%
                dplyr::count(month) %>%
                dplyr::filter(n == max(n)) %>%
                dplyr::pull(month)

            if(length(month_summer) == 0) month_summer <- 1

            month_winter <- site_data %>%
                dplyr::mutate(month = lubridate::month(StartDate)) %>%
                dplyr::filter(Season == "W", !is.na(month)) %>%
                dplyr::count(month) %>%
                dplyr::filter(n == max(n)) %>%
                dplyr::pull(month)

            if(length(month_winter) == 0) month_winter <- 7

            site_data_w_miss <- site_data_w_miss %>%
                dplyr::mutate(StartDate = dplyr::case_when(is.na(StartDate) & Season == "S" ~ as.Date(paste(Year, month_summer, "01", sep = "-")),
                                                           is.na(StartDate) & Season == "W" ~ as.Date(paste(Year, month_winter, "01", sep = "-")),
                                                           TRUE ~ StartDate) ) %>%
                dplyr::arrange(StartDate, Season)

            # Save to temp file
            saveRDS(site_data_w_miss, file = temp_file)

        }

        # Prepare counts for the species and site of interest
        counts[[s]] <- BIRDIE::prepSsmData(counts = site_data_w_miss,
                                           spp_sel = sp_code,
                                           keep = c("LocationCode", "CountCondition", "X", "Y", "Year"))

    }

    # Row-bind all sites
    counts <- dplyr::bind_rows(counts)

    # Save to disk
    #saveRDS(counts, file.path(config$data_outdir, sp_code, paste0(sp_code, "_cwac_data_w_miss.rds")))


    # Annotate with GEE -------------------------------------------------------

    print("Annotating with covariates from GEE")

    #counts <- readRDS(file.path(config$data_outdir, sp_code, paste0(sp_code, "_cwac_data_w_miss.rds")))

    prepGEESpCountData(counts, sp_code, catchment, ...)

}
