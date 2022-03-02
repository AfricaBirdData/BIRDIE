#' Prepare data for fitting a model with occuR
#'
#' @description This function takes generic occupancy visit and site data
#' without species observations, adds detection data for a species in a given
#' year and formats the output to match occuR standards.
#' @param spp_code SAFRING_No of the species of interest as extracted from ABAP.
#' Ignored if download set to FALSE.
#' @param years A numeric vector with elements corresponding to the years we
#' want data for. Ignored if download set to FALSE.
#' @param site_data A dataframe with occupancy site data.
#' @param visit_data A dataframe with occupancy visit data.
#' @param download Indicates whether ABAP data must be downloaded for the species
#' and years indicated by 'spp_code' and 'years'. If TRUE, data will be downloaded
#' from ABAP once per session and cached in a temp file. After this the cached
#' file will be used, unless download is set to "force", in which case data will
#' be downloaded regardless of the cached file. If FALSE data will not be
#' downloaded or retrieved from cache. Note that in this case there won't be
#' any detection data! Only sites and associated covariates.
#' @param save If TRUE data is saved to disc, but see 'overwrite'.
#' @param overwrite A character vector with the data that should be overwriten in
#' case it is already present on disc. It can be any combination of
#' c("site", "visit", "det"), with site referring to site data, visit to visit data
#' and det to detection data. Site and visit data are typically common for
#' multiple species and we might not want to save it all the time.
#' @param config A list of configuration parameters see \link{configPreambOccuR}
#'
#' @return A list containing two data frames: one with site data and one with
#' visit data ready to use with the occuR package
#' @export
#'
#' @examples
prepDataOccuR <- function(spp_code = NULL, years = NULL,
                          site_data, visit_data,
                          download = FALSE, save = TRUE,
                          overwrite = NULL, config){

    if(isFALSE(download)){
        if(!is.null(spp_code) | !is.null(years)){
            warning("spp_code and years are ignored if download is FALSE")
        }
    }

    if(!isFALSE(download)){
        if(!is.null(spp_code) & !is.null(years)){

            # Cache file name
            cachefile <- file.path(tempdir(), paste(c(spp_code, years, ".rds"), collapse = "_"))

            if(file.exists(cachefile) & download != "force"){
                print("Using cached file")
                sp_detect <- readRDS(cachefile)

            } else {

                # Download species detection
                print("Downloading from SABAP")

                sp_detect <- ABAP::getAbapData(.spp_code = spp_code,
                                               .region_type = "country",
                                               .region = "South Africa",
                                               .years = years)

                # Save to cache
                saveRDS(sp_detect, cachefile)
            }

            # Add detections to visit data
            visit_data <- visit_data %>%
                dplyr::left_join(sp_detect %>%
                                     dplyr::select(CardNo, StartDate, Pentad, obs = Spp) %>%
                                     dplyr::mutate(obs = ifelse(obs == "-", 0, 1)),
                                 by = c("CardNo", "StartDate", "Pentad"))
        } else {
            stop("spp_code and years must be provided if download is TRUE")
        }
    }

    if(any(is.na(visit_data$obs))){
        stop("NA found in detection data")
    }

    # Long format for site variables and years
    site_data <- site_data %>%
        BIRDIE::gatherYearFromVars(vars = setdiff(names(.),
                                                  c("Pentad", "lon", "lat", "site", "watocc_ever", "dist_coast")),
                                   sep = "_") %>% #  # HARD CODED. Check that these are the variables that don't change over time
        dplyr::mutate(tdiff = tmmx - tmmn)

    # Visited sites
    visits <- visit_data %>%
        dplyr::distinct(Pentad, year) %>%
        dplyr::mutate(keep = 1)

    sites <- unique(visit_data$Pentad)

    # Keep those sites that appear in visits data and drop geometry
    site_data <- site_data %>%
        dplyr::filter(Pentad %in% sites) %>%
        dplyr::group_by(Pentad) %>%
        dplyr::mutate(site = dplyr::cur_group_id()) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(year) %>%
        dplyr::mutate(occasion = cur_group_id()) %>%
        dplyr::ungroup()

    if(any(is.na(site_data$year))){
        warning("There might be covariates that don't change over time other than watocc_ever and dist_coast.")
    }

    # Keep only sites and occasions in visits
    site_data <- site_data %>%
        dplyr::left_join(visits, by = c("Pentad", "year")) %>%
        dplyr::filter(keep == 1) %>%
        dplyr::select(-keep)

    # Transfer site and occasion indicators from site data
    visit_data <- visit_data %>%
        dplyr::left_join(site_data %>%
                             dplyr::select(Pentad, site) %>%
                             dplyr::distinct(),
                         by = "Pentad") %>%
        dplyr::left_join(site_data %>%
                             dplyr::select(year, occasion) %>%
                             dplyr::distinct(),
                         by = "year")

    # Add visit indicator
    visit_data <- visit_data %>%
        dplyr::group_by(site, occasion) %>%
        dplyr::mutate(visit = row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(CardNo, Date))

    if(any(is.na(visit_data$site))){
        warning("Sites differ between site data and visit data!")
    }

    # Remove data from missing pentads? (THIS SHOULDN'T HAPPEN WHEN PENTADS ARE DOWNLOADED FROM THE API)
    visit_data <- visit_data %>%
        dplyr::filter(!is.na(site))

    site_data <- site_data %>%
        dplyr::filter(site %in% unique(visit_data$site))


    # Save to disc ------------------------------------------------------------

    if(save){

        visitfile <- file.path(config$fit_dir, paste0("occur_visit_dat_sa_", config$years_ch, ".csv"))
        sitefile <- file.path(config$fit_dir, paste0("occur_site_dat_sa_", config$years_ch, ".csv"))
        detfile <- file.path(config$fit_dir, spp_code, paste0("occur_det_dat_sa_", config$years_ch, ".csv"))

        if((!file.exists(visitfile)) | (file.exists(visitfile) & ("visit" %in% overwrite))){
            visit_data %>%
                dplyr::select(-obs) %>%
            write.csv(visitfile, row.names = FALSE)
        } else {
            warning("Visit file not saved because 'visit' not in overwrite")
        }

        if((!file.exists(sitefile)) | (file.exists(sitefile) & ("site" %in% overwrite))){
            site_data %>%
                write.csv(sitefile, row.names = FALSE)
        } else {
            warning("Site file not saved because 'site' not in overwrite")
        }

        if((!file.exists(detfile)) | (file.exists(detfile) & ("det" %in% overwrite))){
            visit_data %>%
                dplyr::select(Pentad, year, visit, obs) %>%
                write.csv(detfile, row.names = FALSE)
        } else {
            warning("Detection file not saved because 'det' not in overwrite")
        }

    }

    return(list(site = site_data,
                visit = visit_data))

}
