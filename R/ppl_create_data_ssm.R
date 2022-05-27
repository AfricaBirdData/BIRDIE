#' Prepare covariate data for abundance pipeline
#'
#' @inheritParams ppl_run_pipe_abu1
#' @param sp_code SAFRING reference number of the species we want to analyze.
#' @param catchment A sf object with polygons corresponding to the catchment, or
#' reference area considered for the covariates associated with CWAC sites.
#' @param steps A character vector expressing the processing steps for the CWAC
#' data. It can be one or more of: "missing" - add missing counts as missing
#' data, "gee" - annotate data with covariates from Google Earth Engine, and
#' "subset" - subset data to those sites with presence of the species and
#' a coverage of at least 10 years from 1993 to 2021. It defaults to all steps.
#' @param ... Other arguments passed on to \link{prepGEESpCountData}
#'
#' @return
#' @export
#'
#' @examples
ppl_create_data_ssm <- function(sp_code, year, catchment, config,
                                steps = c("missing", "gee", "subset"), ...){


    # Load species CWAC data
    if("missing" %in% steps | "subset" %in% steps){

        # Download species data
        print("Downloading from CWAC")

        sp_data <- CWAC::getCwacSppCounts(sp_code)

        # Add DuToit's Doug's extra data
        dutoit <- utils::read.csv(file.path(config$data_dir, "28462448_data_2022_doug.csv"))

        # Transfer column types and bind data frames (the below comes from the CWAC package)
        dutoit <- dutoit %>%
            readr::type_convert(col_types = readr::cols(
                .default = readr::col_integer(),
                LocationName = readr::col_character(),
                Province = readr::col_character(),
                Country = readr::col_character(),
                StartDate = readr::col_date(format = ""),
                Season = readr::col_character(),
                TimeStart = readr::col_time(format = ""),
                TimeEnd = readr::col_time(format = ""),
                WetlandThreat = readr::col_logical(),
                Notes = readr::col_character(),
                record_status = readr::col_character(),
                Survey_notes = readr::col_logical(),
                WetIntCode = readr::col_character(),
                Odr = readr::col_character(),
                Family = readr::col_character(),
                Genus = readr::col_character(),
                Species = readr::col_character(),
                Common_group = readr::col_character(),
                Common_species = readr::col_character(),
                Y = readr::col_double(),
                X = readr::col_double()
            ))

        sp_data <- dplyr::bind_rows(sp_data,
                                    dutoit %>%
                                        dplyr::filter(SppRef == sp_code))

        # Clean any data from before 1993 (those should be errors)
        sp_data <- sp_data %>%
            dplyr::filter(Year > 1992)

    }


    # Add missing counts ------------------------------------------------------

    if("missing" %in% steps){

        # In what South African sites was the species present?
        sp_sites <- sp_data %>%
            dplyr::filter(Country == "South Africa",
                          Year %in% config$years) %>%
            dplyr::pull(LocationCode) %>%
            unique()


        # Add missing counts ------------------------------------------------------

        print("Adding missing counts to all sites. This might take a while...")

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

                # Deal with DuToit's extra data
                if(site_sel == 28462448){
                    site_data <- dplyr::bind_rows(site_data,
                                                  dutoit)
                }

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

        attr(counts, "missing") <- TRUE

        # Save to disk
        saveRDS(counts, file.path(tempdir(), paste0(sp_code, "_", config$years_ch, "_cwac_data_w_miss.rds")))
    }


    # Annotate with GEE -------------------------------------------------------

    if("gee" %in% steps){

        outfile <- file.path(tempdir(), paste0(sp_code, "_", config$years_ch, "_cwac_data_w_miss.rds"))

        if(!exists("counts") & file.exists(outfile)){
            counts <- readRDS(outfile)
        } else if(!exists("counts") & !file.exists(outfile)){
            stop("No dataset with missing counts found. Perhaps you need to run the 'missing' step?")
        }

        print("Annotating with covariates from GEE")

        counts <- prepGEESpCountData(counts, sp_code, catchment, config, ...)

        attr(counts, "gee") <- TRUE

    }


    # Subset sites ------------------------------------------------------------

    if("subset" %in% steps){

        if(!exists("counts") || !isTRUE(attr(counts, "gee"))){
            counts <- utils::read.csv(file.path(config$out_dir, sp_code, paste0("abu_gee_data_", sp_code, "_", config$years_ch, ".csv")))
        }

        print("Finding suitable CWAC sites (>= 15 year coverage from 1993 to 2021)")

        # Find those sites that have a coverage of at
        sites_good <- sp_data %>%
            dplyr::count(LocationCode) %>%
            dplyr::filter(n > 10) %>%
            dplyr::pull(LocationCode)

        counts <- counts %>%
            dplyr::filter(LocationCode %in% sites_good)

        outfile <- file.path(config$out_dir, sp_code, paste0("abu_model_data_", sp_code, "_", config$years_ch, ".csv"))

        utils::write.csv(counts, outfile, row.names = FALSE)

        print(paste("Final counts dataset saved at", outfile))

    }

    return(counts)

}