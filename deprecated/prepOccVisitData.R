#' Prepare occupancy visit data
#'
#' @param region A character string defining the region of interest. At the
#' moment only countries are supported. SABAP pentads contained in the region
#' will be used as sites.
#' @param sites An sf object with the sites we are interested in. If provided,
#' 'region' is ignored.
#' @param species The SABAP numeric code of the species of interest.
#' @param years A numeric vector of years we want to estimate occupancy for.
#' Covariates for these years will be appended to the data.
#' @param clim_covts A character vector with the abbreviation of the climatic
#' covariates we want to append to the data. At the moment, the Terraclimate
#' data set is used. See \link[climateR]{getTerraClim} for the variables
#' available.
#' @param covts_dir Directory where covariate data are found.
#' @param file_fix A character vector with two elements corresponding to
#' additional characters that are found before and after clim_covts.
#' @param downld_path If provided administrative boundaries layer will be looked
#' for and stored at here.
#' @param savedir Directory where sites should be saved at before covariate
#' extraction. Intersecting region with SABAP pentads is time consuming, so you
#' may want to save the selected pentads as soon as possible. The directory
#' must include the name and extension of the file (.rds).
#'
#' @return
#' @export
#'
#' @examples
prepOccVisitData <- function(region = NULL, sites, species, years,
                             clim_covts, covts_dir, file_fix, downld_path = NULL,
                             savedir = NULL){

    if("sequential:" %in% utils::capture.output(future::plan())){
        print("Running on a single core. Set future::plan('multisession') to run on multiple cores")
    }

    # Prepare occupancy visit data -------------------------------------------

    print("Downloading SABAP data. This might take a while...")

    # Extract SABAP data
    if(is.null(region)){
        pa_dat <- SABAP::getSabapData(.spp_code = species,
                                      .region_type = "pentad",
                                      .region = unique(sites$Name),
                                      .years = years)
    } else {
        pa_dat <- SABAP::getSabapData(.spp_code = species,
                                      .region_type = "country",
                                      .region = region,
                                      .years = years)
    }

    # Add year and month
    pa_dat <- pa_dat %>%
        dplyr::mutate(year = lubridate::year(StartDate),
                      month = lubridate::month(StartDate))

    # Filter out those pentads that are not in the region
    pa_dat <- pa_dat %>%
        dplyr::filter(Pentad %in% unique(sites$Name))

    # Join detection data and geographic data
    pa_dat <- pa_dat %>%
        dplyr::left_join(sites %>%
                             sf::st_drop_geometry() %>%
                             dplyr::select(Name, lon, lat, id),
                         by = c("Pentad" = "Name"))

    # Fix names
    pa_dat <- pa_dat %>%
        dplyr::mutate(obs = dplyr::if_else(Spp == "-", 0, 1)) %>%
        dplyr::select(id, CardNo, lon, lat, year, month, StartDate, EndDate,
                      Pentad, TotalHours, obs)


    # Annotate occupancy visit data with climatic covariates

    print("Extracting climatic covariates")

    for(i in seq_along(clim_covts)){
        pa_dat <- BIRDIE::addOccVisitCovt(visits = pa_dat,
                                          sites = sites,
                                          covt = clim_covts[i],
                                          covts_dir = covts_dir,
                                          file_fix = file_fix)
    }


    # Prepare non-linear terms ------------------------------------------------

    # # Load detection/non-detection data
    # pa_dat <- readRDS(paste0("analysis/data/visit_dat_", sp_sel, "_wcovts.rds"))
    #
    # # Add cyclic basis values for time of year to visit data
    # spl_bs <- mgcv::cSplineDes(x = pa_dat$month,
    #                            knots = seq(1, 12, length.out = 5),
    #                            ord = 4, derivs = 0)
    # colnames(spl_bs) <- paste0("toy.", 1:ncol(spl_bs))
    # pa_dat <- cbind(pa_dat, spl_bs)
    #
    # saveRDS(pa_dat, paste0("analysis/data/visit_dat_", sp_sel, "_wcovts.rds"))

    return(pa_dat)

}
