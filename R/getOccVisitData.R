getOccVisitData <- function(region_type, region, species, ...){

    # Find pentads in the region
    pentads_sel <- readRDS("analysis/data/pentads_nw.rds")

    # OR

    if(region_type == "province"){
        country <- "South Africa"
        province <- region
    }

    sf::sf_use_s2(FALSE) # s2 intersection takes very long
    pentads_sel <- getRegionPentads(.country = country, .province = province, ...)

    # Add centroid coordinates
    sf::sf_use_s2(TRUE)
    pentads_sel <- pentads_sel %>%
        dplyr::mutate(lon = sf::st_coordinates(sf::st_centroid(.))[,1],
                      lat = sf::st_coordinates(sf::st_centroid(.))[,2])

    attr(pentads_sel$lon, "names") <- NULL
    attr(pentads_sel$lat, "names") <- NULL

    # Correct potential name problems and get SABAP data
    # region <- tolower(region)
    # region <- gsub(" ", "", region)

    pa_dat <- SABAP::getSabapData(.spp_code = species,
                                  .region_type = region_type,
                                  .region = region)

    # Add year and month
    pa_dat <- pa_dat %>%
        dplyr::mutate(year = lubridate::year(StartDate),
               month = lubridate::month(StartDate))

    # Filter out those pentads that are not in the region
    pa_dat <- pa_dat %>%
        dplyr::filter(Pentad %in% unique(pentads_sel$Name))

    # Join detection data and geographic data
    pa_dat <- pa_dat %>%
        dplyr::left_join(pentads_sel %>%
                             sf::st_drop_geometry() %>%
                             dplyr::select(Name, lon, lat),
                         by = c("Pentad" = "Name"))

    # Add id and fix names. Names follow Rcppocc requirements
    pa_dat <- pa_dat %>%
        dplyr::mutate(id = dplyr::row_number(),
                      PAdata = dplyr::if_else(Spp == "-", 0, 1)) %>%
        dplyr::select(id, CardNo, lon, lat, year, month, StartDate, EndDate,
                      SiteName = Pentad, TotalHours, PAdata)

    # There are no covariates for 2020 so...
    pa_dat <- pa_dat %>%
        dplyr::filter(year < 2020)

    return(list(visits = pa_dat, sites = pentads_sel))

}
