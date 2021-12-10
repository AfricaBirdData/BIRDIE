#' Prepare occupancy site data
#'
#' @param region A character string defining the region of interest. At the
#' moment only countries are supported. SABAP pentads contained in the region
#' will be used as sites.
#' @param sites An sf object with the sites we are interested in. If provided,
#' 'region' is ignored.
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
prepOccSiteData <- function(region = NULL, sites = NULL, years, clim_covts,
                            covts_dir, file_fix, downld_path = NULL,
                            savedir = NULL){

    if("sequential:" %in% utils::capture.output(future::plan())){
        print("Running on a single core. Set future::plan('multisession') to run on multiple cores")
    }

    # Extract region pentads if necessary ------------------------------------

    if(is.null(sites)){

        print("Intersecting SABAP pentads. This might take some time.")

        sf::sf_use_s2(FALSE) # s2 intersection takes very long

        sites <- SABAP::getRegionPentads(.country = region, .path = downld_path)

        if(!is.null(savedir)){

            print(paste0("Saving sites at ", savedir, "sites.rds"))

            saveRDS(sites, paste0(savedir, "sites.rds"))

        }

    }

    # Add centroids -----------------------------------------------------------

    # Add centroid coordinates

    print("Calculating centroids")

    sf::sf_use_s2(TRUE)
    cc <- sites %>%
        sf::st_centroid() %>%
        sf::st_coordinates()

    attr(cc[,1], "names") <- NULL
    attr(cc[,2], "names") <- NULL

    sites <- sites %>%
        dplyr::mutate(lon = cc[,1],
                      lat = cc[,2]) %>%
        dplyr::arrange(lat, lon) %>%
        dplyr::mutate(id = dplyr::row_number())


    # Prepare site covariate data --------------------------------------------

    print("Extracting climatic covariates")

    for(i in seq_along(clim_covts)){
        sites <- BIRDIE::addOccSiteCovt(sites,
                                        covt = clim_covts[i],
                                        years = years,
                                        covts_dir = covts_dir,
                                        file_fix = file_fix)
    }


    # Prepare surface water ---------------------------------------------------

    print("Extracting surface water")

    if(length(grep("surf_water", list.files(covts_dir))) == 0){
        stop("Please name surface water layers 'surf_water_xxE_xxS.tif")
    }

    # Define coordinates of surface water maps
    cc <- list(c(10, 20),
               c(20, 20),
               c(30, 20),
               c(10, 30),
               c(20, 30),
               c(30, 30))

    out <- vector("list", length = length(cc))

    for(i in seq_along(cc)){

        # Load water raster
        water <- raster::raster(paste0(covts_dir, "surf_water_", cc[[i]][1], "E_", cc[[i]][2], "S.tif"))

        water <- raster::crop(water, sites)

        # Reclassify 255 to NA
        water <- raster::reclassify(water, rcl = matrix(c(255, NA), ncol = 2))

        # Extract values
        sites_water <- BIRDIE::exactExtractParll(water, sites,
                                                 ncores = future::nbrOfWorkers(),
                                                 fun = "mean")

        sites_water <- sites_water %>%
            dplyr::rename(water = vals)

        out[[i]] <- sites_water %>%
            dplyr::filter(!is.na(water))

        sites <- sites_water %>%
            dplyr::filter(is.na(water)) %>%
            dplyr::select(-water)

    }

    sites <- do.call("rbind", out)
    sites <- dplyr::arrange(sites, id)

    return(sites)

}
