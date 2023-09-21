#' Prepare covariate data for abundance pipeline
#'
#' @inheritParams ppl_run_pipe_abu1
#' @param sp_code SAFRING reference number of the species we want to analyse.
#' @param year Year for which the SSM data should be prepared.
#' @param catchment A sf object with polygons corresponding to the catchment, or
#' reference area considered for the covariates associated with CWAC sites.
#' @param steps A character vector expressing the processing steps for the CWAC
#' data. It can be one or more of: "missing" - add missing counts as missing
#' data, "gee" - annotate data with covariates from Google Earth Engine,
#' "subset" - subset data to those sites with presence of the species and
#' a coverage of at least 10 years from 1993 to 2021, and "model" prepare data
#' for model fitting by adding some aux variables and ordering the data.
#' It defaults to all steps.
#' @param force_gee Logical. If TRUE (default), then count data will be annotated with
#' environmental information from GEE. If FALSE, then we assume data have already been
#' annotated and only some manipulation to join catchment environmental data with
#' count data is required.
#' @param force_gee_upload Logical. If TRUE (default), then the catchment polygons
#' will be uploaded to GEE under the name 'quin_catchm'. If FALSE, then we assume
#' these polygons are already #' present in GEE server and we don't need to upload them again.
#' @param ... Other arguments passed on to \link{prepGEECatchmData}
#'
#' @return
#' @export
#'
#' @examples
ppl_create_data_ssm <- function(sp_code, year, catchment, config, force_gee = TRUE,
                                force_gee_upload = TRUE,
                                steps = c("subset", "missing", "gee", "model"),
                                ...){

    varargs <- list(...)

    # Load species CWAC data
    if("missing" %in% steps | "subset" %in% steps){

        # Download species data
        message("Downloading from CWAC")

        sp_data <- CWAC::getCwacSppCounts(sp_code)

        # Add DuToits data
        sp_data <- addDuToitCounts(sp_data, config)

        # Clean any data from before 1993 (those should be errors)
        sp_data <- sp_data %>%
            dplyr::filter(Year > 1992,
                          SppRef == sp_code)

        # Subset data to the country of interest
        region <- dplyr::case_when(config$region == "ZA" ~ "South Africa",
                                   config$region == "KE" ~ "Kenya")

        sp_data <- sp_data %>%
            dplyr::filter(Country == region)

        if(nrow(sp_data) == 0){
            warning(paste("There seems to be no data for species", sp_code))
            sink(setSpOutFilePath("No_CWAC_data", config, config$years_ch, sp_code, ".txt"))
            message(paste("No CWAC data for species", sp_code, "on years 1993-2021"))
            sink()
            return(1)
        }

    }


    # Subset sites ------------------------------------------------------------

    if("subset" %in% steps){

        message("Finding suitable CWAC sites (>= 10 year coverage from 1993 to 2021)")

        # Find those sites where the species was detected in at least 10 years between 1993-2021
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
        if(length(sites_good) < 1){

            sink(setSpOutFilePath("No_CWAC_sites", config, config$years_ch, sp_code, ".txt"))
            message(paste("No suitable CWAC sites for species", sp_code, "on years 1993-2021"))
            sink()

            # return(1)

        } else if(length(sites_good) < 5){

            sink(setSpOutFilePath("Less_5_CWAC_sites", config, config$years_ch, sp_code, ".txt"))
            message(paste("Less than 5 suitable CWAC sites for species", sp_code, "on years 1993-2021"))
            sink()

        }

        sp_data_sel <- sp_data %>%
            dplyr::mutate(good_site = ifelse(LocationCode %in% sites_good, 1, 0))

        # Save to disk at temporary location
        tmp_sp_data_sel <- file.path(tempdir(), paste0(sp_code, "_", config$years_ch, "_", config$region, "_cwac_data_subset.rds"))
        saveRDS(sp_data_sel, tmp_sp_data_sel)

        # Save to disk permanently?
        # saveRDS(counts, file.path("analysis/out_nosync", sp_code, paste0(sp_code, "_", config$years_ch, "_cwac_data_subset.rds")))

    }



    # Add missing counts ------------------------------------------------------

    if("missing" %in% steps){

        message("Adding missing counts to all sites. This might take a while...")

        tmp_sp_data_sel <- file.path(tempdir(), paste0(sp_code, "_", config$years_ch, "_", config$region, "_cwac_data_subset.rds"))

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
            dplyr::filter(Year %in% config$years) %>%
            dplyr::pull(LocationCode) %>%
            unique()

        # List to store data from multiple sites
        counts <- vector("list", length = length(sp_sites))

        # For each site:
        for(s in seq_along(sp_sites)){

            site_sel <- sp_sites[s]

            good <- sp_data_sel %>%
                dplyr::filter(LocationCode == site_sel) %>%
                dplyr::pull(good_site) %>%
                unique()

            # Temp file to save to/retrieve from. We use this temporary file because we
            # are not concerned with any particular species when including missing
            # counts. And because adding missing counts is a time-consuming process,
            # we can make use of a pre-computed temp file that can be used multiple times.
            temp_file <- file.path(tempdir(), paste(site_sel, config$years_ch, "cwac_dat.rds", sep = "_"))

            if(file.exists(temp_file)){

                site_data_w_miss <- readRDS(temp_file)

            } else {

                # Download site data and prepare missing visits
                site_data <- CWAC::getCwacSiteCounts(site_sel)

                # Deal with DuToit's extra data
                if(site_sel == 28462448){
                    site_data <- addDuToitCounts(site_data, config)
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
                                                               TRUE ~ StartDate)) %>%
                    dplyr::arrange(StartDate, Season)

                # Save to temp file
                saveRDS(site_data_w_miss, file = temp_file)

            }

            # Add information as to whether this site is good for the species. Note that
            # this has to be outside of the if statement, because it needs to be evaluated
            # for each species. So it cannot use the cached version (temp_file).
            site_data_w_miss <- site_data_w_miss %>%
                dplyr::mutate(good_site = good)

            # Prepare counts for the species and site of interest
            counts[[s]] <- prepSsmData(counts = site_data_w_miss,
                                       spp_sel = sp_code,
                                       keep = c("Country", "LocationCode", "CountCondition", "X", "Y", "Year", "good_site"))

        }

        # Row-bind all sites
        counts <- dplyr::bind_rows(counts)

        attr(counts, "missing") <- TRUE

        # Save to disk at temporary location
        utils::write.csv(counts,
                         setSpOutFilePath("cwac_data_w_miss", config, config$years_ch, sp_code, ".csv"),
                         row.names = FALSE)

        # Save to disk permanently?
        # saveRDS(counts, file.path("analysis/out_nosync", sp_code, paste0(sp_code, "_", config$years_ch, "_cwac_data_w_miss.rds")))
    }


    # Annotate with GEE -------------------------------------------------------

    if("gee" %in% steps){

        outfile <- setSpOutFilePath("cwac_data_w_miss", config, config$years_ch, sp_code, ".csv")

        if(!exists("counts") & file.exists(outfile)){
            counts <- utils::read.csv(outfile, colClasses = c(LocationCode = "character"))
        } else if(!exists("counts") & !file.exists(outfile)){
            stop("No dataset with missing counts found. Perhaps you need to run the 'missing' step?")
        }

        # We annotate only good sites
        counts <- counts %>%
            dplyr::filter(good_site == 1) %>%
            dplyr::mutate(id_count = dplyr::row_number())

        # There are catchment data only for South Africa ATM, so they are the only GEE annotated data
        if(config$region == "ZA" && nrow(counts) > 0){

            message("Annotating with covariates from GEE")

            geefile <- file.path(config$out_dir, paste0("catchm_dat_sa_gee_", config$years_ch, ".csv"))

            if(!force_gee & file.exists(geefile)){
                warning("File with covariates already on disk. Set force_gee = TRUE to overwrite")
            }

            if(!file.exists(geefile) | force_gee){

                rgee::ee_check()
                rgee::ee_Initialize()

                if(force_gee_upload){
                    catchment %>%
                        ABDtools::uploadFeaturesToEE(asset_id = file.path(rgee::ee_get_assethome(), 'quin_catchm'),
                                                     load = FALSE,
                                                     monitor = varargs$monitor)
                }

                gee_catchm <- prepGEECatchmData(sp_code, catchment, config, monitor = varargs$monitor)

                outfile <- file.path(config$out_dir, paste0("catchm_dat_sa_gee_", config$years_ch, ".csv"))

                utils::write.csv(gee_catchm, outfile, row.names = FALSE)

                message(paste("Catchment data with GEE covts saved at", outfile))

            } else {
                gee_catchm <- utils::read.csv(geefile)
            }

            # Find nearest catchment to site counts. Nearest because some sites are on the
            # coast potentially falling offshore
            counts <- counts %>%
                sf::st_as_sf(coords = c("X", "Y"), dim = "XY", crs = sf::st_crs(4326))

            int_index <- counts %>%
                sf::st_nearest_feature(catchment)

            counts <- counts %>%
                dplyr::mutate(UNIT_ID = catchment$UNIT_ID[int_index])

            counts <- counts %>%
                sf::st_drop_geometry()

            names_long <- setdiff(names(gee_catchm)[grep("_", names(gee_catchm))], c("UNIT_ID", "watocc_ever", "dist_coast"))

            catchm_long <- gatherYearFromVars(gee_catchm, names_long, sep = "_")

            # Before joining we add one year to environmental data. This is because
            # summer waterbird populations should be affected by the conditions in the previous
            # year rather than by conditions in the following year.
            catchm_long$year <- catchm_long$year + 1

            counts <- counts %>%
                dplyr::left_join(catchm_long, by = c("UNIT_ID", "year"))

            attr(counts, "gee") <- TRUE

        }

        if(nrow(counts) > 0){

            # Set species code and output file name
            outfile <- setSpOutFilePath("abu_gee_data", config, config$years_ch, sp_code, ".csv")

            counts %>%
                utils::write.csv(outfile, row.names = FALSE)

            message(paste("Dataset with available GEE covts saved at", outfile))

        }

    }


    # Prepare data for modelling ----------------------------------------------

    if("model" %in% steps){

        outfile <- setSpOutFilePath("abu_gee_data", config, config$years_ch, sp_code, ".csv")

        if(file.exists(outfile)){

            counts <- read.csv(outfile, colClasses = c(LocationCode = "character"))

            # Remove seasons other than summer and winter
            counts_mod <- counts %>%
                dplyr::filter(Season %in% c("S", "W"))

            # Create other useful variables
            counts_mod <- counts_mod %>%
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
                dplyr::arrange(site_id, year_id, season_id, visit_id) %>%
                dplyr::mutate(spp_id = sp_code)

            # Order data by: site, year, season
            counts_mod <- counts_mod %>%
                dplyr::arrange(site_id, year_id, season_id)

            return(counts_mod)

        } else {
            warning("No data suitable for modelling found")
            return(1)
        }

    }

}
