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
#' downloaded or retrieved from cache.
#'
#' @return A list containing two data frames: one with site data and one with
#' visit data ready to use with the occuR package
#' @export
#'
#' @examples
prepDataOccuR <- function(spp_code = NULL, years = NULL,
                          site_data, visit_data,
                          download = FALSE){

    if(isFALSE(download)){
        if(!is.null(spp_code) | !is.null(years)){
            warning("spp_code and years are ignored if download is FALSE")
        }
    }

    if(!isFALSE(download)){
        if(!is.null(spp_code) & !is.null(years)){

            # Cache file name
            cachefile <- tempfile(pattern = paste(c("abap", spp_code, years), collapse = "_"),
                                  fileext = ".rds")

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
                                     dplyr::select(CardNo, obs = Spp) %>%
                                     dplyr::mutate(obs = if_else(obs == "-", 0, 1)),
                                 by = "CardNo")
        } else {
            stop("spp_code and years must be provided if download is TRUE")
        }
    }

    # Long format for site variables and years
    site_data <- site_data %>%
        sf::st_drop_geometry() %>%
        BIRDIE::gatherYearFromVars(vars = names(.)[-c(1:5)], sep = "_") %>% #  # HARD CODED. Check that 1:5 are the variables that don't change over time
        dplyr::mutate(tdiff = tmmx - tmmn)

    # Visited sites
    visits <- visit_data %>%
        dplyr::distinct(Pentad, year) %>%
        dplyr::mutate(keep = 1)

    sites <- unique(visit_data$Pentad)

    # Keep those sites that appear in visits data and drop geometry
    site_data <- site_data %>%
        dplyr::filter(Name %in% sites) %>%
        dplyr::group_by(Name) %>%
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
        dplyr::left_join(visits, by = c("Name" = "Pentad", "year")) %>%
        dplyr::filter(keep == 1) %>%
        dplyr::select(-keep)

    # Transfer site and occasion indicators from site data
    visit_data <- visit_data %>%
        dplyr::left_join(site_data %>%
                             dplyr::select(Name, site) %>%
                             dplyr::distinct(),
                         by = c("Pentad" = "Name")) %>%
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

    # Make data.tables
    site_data <- data.table::as.data.table(site_data)
    visit_data <- data.table::as.data.table(visit_data)


    # Scaling -----------------------------------------------------------------

    scaling <- list(visit = NULL,
                   site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")) # HARD CODED

    if(!is.null(scaling)){

        if(!is.null(scaling$visit)){
            visit_data <- visit_data %>%
                dplyr::mutate(dplyr::across(.col = dplyr::all_of(scaling$visit), .fns = ~scale(.x)))
            attr(visit_data, "scaling") <- TRUE
        } else {
            attr(visit_data, "scaling") <- FALSE
        }

        if(!is.null(scaling$site)){
            site_data <- site_data %>%
                dplyr::mutate(dplyr::across(.col = dplyr::all_of(scaling$site), .fns = ~scale(.x)))
            attr(site_data, "scaling") <- TRUE
        } else {
            attr(site_data, "scaling") <- FALSE
        }

    }

    # Remove data from missing pentads? (THIS SHOULDN'T HAPPEN WHEN PENTADS ARE DOWNLOADED FROM THE API)
    visit_data <- visit_data %>%
        dplyr::filter(!is.na(site))

    site_data <- site_data %>%
        dplyr::filter(site %in% unique(visit_data$site))

    return(list(site = site_data,
                visit = visit_data))

}
