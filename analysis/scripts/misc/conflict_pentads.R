# Identify conflicting pentads --------------------------------------------

# Any pentads in occu_nw that not in pentads_nw?
any(!unique(occu_nw$Pentad) %in% unique(pentads_nw$Name)) # YES, WHY? Check with Michael

# Which pentands in occu_nw not in pentads_nw? Viceversa?
sab_not_geo <- unique(occu_nw$Pentad)[!unique(occu_nw$Pentad) %in% unique(pentads_nw$Name)]
geo_not_sab <- unique(pentads_nw$Name)[!unique(pentads_nw$Name) %in% unique(occu_nw$Pentad)]

# Plot
pentads %>%
    filter(Name %in% occu_nw$Pentad) %>%
    st_geometry() %>%
    plot()

pentads %>%
    filter(Name %in% pentads_nw$Name) %>%
    st_geometry() %>%
    plot(col = "green", lwd = 1, add = T)

# Add region
region <- raster::getData("GADM", download = TRUE, country = "South Africa",
                          level = 1, path = "analysis/data") %>%
    st_as_sf() %>%
    filter(NAME_1 == "North West")

plot(st_geometry(region), add = T)

ggplot() +
    geom_sf(data = filter(pentads, Name %in% occu_nw$Pentad), fill = "yellow", alpha = 1) +
    geom_sf(data = filter(pentads, Name %in% pentads_nw$Name), fill = "blue", alpha = 0.2) +
    geom_sf(data = region, fill = NA)

