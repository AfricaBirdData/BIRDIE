#' Prepare covariate data for abundance pipeline
#'
#' @inheritParams ppl_run_pipe_abu1
#' @param sp_code SAFRING reference number of the species we want to analyze.
#' @param catchment A sf object with polygons corresponding to the catchment, or
#' reference area considered for the covariates associated with CWAC sites.
#' @param steps A character vector expressing the processing steps for the CWAC
#' data. It can be one or more of: "missing" - add missing counts as missing
#' data, "gee" - annotate data with covariates from Google Earth Engine,
#' "subset" - subset data to those sites with presence of the species and
#' a coverage of at least 10 years from 1993 to 2021, and "model" prepare data
#' for model fitting by adding some aux variables and ordering the data.
#' It defaults to all steps.
#' @param ... Other arguments passed on to \link{prepGEESpCountData}
#'
#' @return
#' @export
#'
#' @examples
ppl_create_data_ssm <- function(sp_code, year, catchment, config,
                                steps = c("missing", "gee", "subset", "model"),
                                ...){


    # Load species CWAC data
    if("missing" %in% steps | "subset" %in% steps){

        # Download species data
        message("Downloading from CWAC")

        sp_data <- CWAC::getCwacSppCounts(sp_code)

        sp_data <- sp_data %>%
            dplyr::select(LocationCode, LocationName, Province, Country, Year,
                          StartDate, Season, SppRef, WetIntCode, Species,
                          Common_group, Common_species, Count, X, Y)

        # Add DuToit's Doug's extra data
        dutoit <- utils::read.csv(file.path(config$data_dir, "28462448_data_2022_doug.csv")) %>%
            dplyr::select(LocationCode, LocationName, Province, Country, Year,
                          StartDate, Season, SppRef, WetIntCode, Species,
                          Common_group, Common_species, Count, X, Y)

        # Transfer column types and bind data frames (the below comes from the CWAC package)
        dutoit <- dutoit %>%
            readr::type_convert(col_types = readr::cols(
                .default = readr::col_integer(),
                LocationName = readr::col_character(),
                Province = readr::col_character(),
                Country = readr::col_character(),
                StartDate = readr::col_date(format = ""),
                Season = readr::col_character(),
                # TimeStart = readr::col_time(format = ""),
                # TimeEnd = readr::col_time(format = ""),
                WetlandThreat = readr::col_logical(),
                # Notes = readr::col_character(),
                # record_status = readr::col_character(),
                # Survey_notes = readr::col_logical(),
                WetIntCode = readr::col_character(),
                # Odr = readr::col_character(),
                # Family = readr::col_character(),
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


    # Subset sites ------------------------------------------------------------

    if("subset" %in% steps){

        message("Finding suitable CWAC sites (>= 10 year coverage from 1993 to 2021)")

        # Find those sites that have a coverage of at least 10 years between 1993-2021
        # Sites must have counts in both summer and winter. And both summer and winter
        # must have been counted at least 5 times each
        sites_good <- sp_data %>%
            # Counted at least 5 times in summer or 5 times in winter
            dplyr::filter(Season %in% c("S", "W"),
                          !is.na(Count),
                          Year %in% 1993:2021) %>%
            dplyr::count(LocationCode, Year, Season) %>%
            dplyr::count(LocationCode, Season) %>%
            dplyr::filter(n > 4) %>%
            # Counted in both seasons (if counted in both and at least 5 in each, then at least 10 in total as well)
            dplyr::count(LocationCode) %>%
            dplyr::filter(n > 1) %>%
            dplyr::pull(LocationCode)

        # If less than 5 sites create notification, because we might need a different model
        if(length(sites_good) < 5){

            sink(setSpOutFilePath("Less_5_sites", config, sp_code, ".txt"))
            message(paste("Less than 5 sites for species", sp_code, "on year", config$year))
            sink()

            return(1)

        }

        sp_data_sel <- sp_data %>%
            dplyr::filter(LocationCode %in% sites_good)

        # Save to disk at temporary location
        tmp_sp_data_sel <- file.path(tempdir(), paste0(sp_code, "_", config$years_ch, "_cwac_data_subset.rds"))
        saveRDS(sp_data_sel, tmp_sp_data_sel)

        # Save to disk permanently?
        # saveRDS(counts, file.path("analysis/out_nosync", sp_code, paste0(sp_code, "_", config$years_ch, "_cwac_data_subset.rds")))

    }



    # Add missing counts ------------------------------------------------------

    if("missing" %in% steps){

        message("Adding missing counts to all sites. This might take a while...")

        tmp_sp_data_sel <- file.path(tempdir(), paste0(sp_code, "_", config$years_ch, "_cwac_data_subset.rds"))

        if(!exists("sp_data_sel")){

            if(file.exists(tmp_sp_data_sel)){
                sp_data_sel <- readRDS(tmp_sp_data_sel)
                message("Using cached subset of suitable sites")
            } else {
                message("CAUTION: adding missing counts to species data without subsetting suitable sites!")
                sp_data_sel <- sp_data
            }

        }

        # In what South African sites was the species present?
        sp_sites <- sp_data_sel %>%
            dplyr::filter(Country == "South Africa",
                          Year %in% config$years) %>%
            dplyr::pull(LocationCode) %>%
            unique()


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

        # Save to disk at temporary location
        saveRDS(counts, file.path(tempdir(), paste0(sp_code, "_", config$years_ch, "_cwac_data_w_miss.rds")))

        # Save to disk permanently?
        # saveRDS(counts, file.path("analysis/out_nosync", sp_code, paste0(sp_code, "_", config$years_ch, "_cwac_data_w_miss.rds")))
    }


    # Annotate with GEE -------------------------------------------------------

    if("gee" %in% steps){

        outfile <- file.path(tempdir(), paste0(sp_code, "_", config$years_ch, "_cwac_data_w_miss.rds"))

        if(!exists("counts") & file.exists(outfile)){
            counts <- readRDS(outfile)
        } else if(!exists("counts") & !file.exists(outfile)){
            stop("No dataset with missing counts found. Perhaps you need to run the 'missing' step?")
        }

        message("Annotating with covariates from GEE")

        geefile <- file.path(config$out_dir, paste0("catchm_dat_sa_gee_", config$years_ch, ".csv"))

        if(!file.exists(geefile) | force_gee){
            prepGEECatchmData(sp_code, catchment, config, ...)
        } else {
            gee_catchm <- utils::read.csv(geefile)
        }

        try(
            counts <- prepGEESpCountData(counts, sp_code, catchment, config, ...)
        )

        attr(counts, "gee") <- TRUE

    }



    # Prepare data for modelling ----------------------------------------------


    if("model" %in% steps){

        outfile <- setSpOutFilePath("abu_gee_data", config, sp_code, ".csv")

        if(!exists("counts") & file.exists(outfile)){
            counts <- read.csv(outfile)
        } else if(!exists("counts") & !file.exists(outfile)){
            stop("No dataset with subset counts found. Perhaps you need to run the 'subset' step?")
        }

        counts_mod <- counts %>%
            dplyr::mutate(season_id = dplyr::case_when(Season == "S" ~ 1,
                                                       Season == "W" ~ 2,
                                                       TRUE ~ 3),
                          date = lubridate::date(StartDate)) %>%
            dplyr::group_by(LocationCode) %>%
            dplyr::mutate(site_id = dplyr::cur_group_id()) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(year) %>%
            dplyr::mutate(year_id = dplyr::cur_group_id()) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(site_id, year_id, season_id) %>%
            dplyr::arrange(date) %>%
            dplyr::mutate(visit_id = dplyr::row_number()) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(site_id, year_id, season_id, visit_id)

        # Remove seasons other than summer and winter
        counts_mod <- counts_mod %>%
            dplyr::filter(season_id != 3)

        # Add counts that come from the same site, season and year
        gen_vars <- counts_mod %>%
            dplyr::select(year, Season, spp, LocationCode, UNIT_ID,
                          season_id, site_id, year_id) %>%
            dplyr::distinct()

        counts_mod <- counts_mod %>%
            dplyr::group_by(site_id, year_id, season_id) %>%
            dplyr::summarise(count = sum(count, na.rm = TRUE),
                             dplyr::across(.cols = c(prcp_mean, tmmn_mean, tmmx_mean, pdsi_mean, watext_count, watrec_mean),
                                           ~mean(.x))) %>%
            dplyr::ungroup()

        counts_mod <- counts_mod %>%
            dplyr::left_join(gen_vars, by = c("site_id", "year_id", "season_id"))


        # Create variables that are change in covariates
        # counts_mod <- counts_mod %>%
        #     dplyr::group_by(site_id, year_id, season_id) %>%
        #     dplyr::mutate(dplyr::across(.cols = c(pdsi_mean, watext_count, watrec_mean),
        #                                 .fns = mean, .names = "mean_{.col}")) %>%
        #     dplyr::group_by(site_id, season_id) %>%
        #     dplyr::mutate(dplyr::across(.cols = c(mean_pdsi_mean, mean_watext_count, mean_watrec_mean),
        #                                 .fns = ~c(diff(.x), NA), .names = "diff_{.col}")) %>%
        #     dplyr::rename_with(~gsub("_mean_", "_", .x), .cols = dplyr::everything()) %>%
        #     dplyr::ungroup()

        # counts_mod <- counts_mod %>%
        #     dplyr::select(year, Season, StartDate, LocationCode, year_id, season_id, site_id, visit_id, id_count) %>%
        #     dplyr::arrange(site_id, year_id, season_id, visit_id)

    }

    return(counts)

}
