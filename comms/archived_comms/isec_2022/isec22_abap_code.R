library(ABAP)
library(sf)
library(ggplot2)
library(dplyr)
library(rgee)

ee_check()
ee_Initialize(drive = TRUE)

flamingos <- searchAbapSpecies("Flamingo")

gf_abap <- getAbapData(86, "province", "Western Cape", 2018)

pentads <- ABAP::getRegionPentads("province", "Western Cape")

plot(st_geometry(pentads))


pentads %>%
    ggplot() +
    geom_sf() +
    theme_void()

ggsave("comms/isec_2022/isec22_abap/pentads_wc.png")


pentads %>%
    left_join(gf_abap %>%
                  mutate(obs = if_else(Spp == "-", 0, 1)) %>%
                  group_by(Pentad) %>%
                  summarize(det = max(obs)) %>%
                  ungroup(),
              by = c("pentad" = "Pentad")) %>%
    ggplot() +
    geom_sf(aes(fill = factor(det))) +
    scale_fill_viridis_d(name = "detection", option = "cividis") +
    theme_void()

ggsave("comms/isec_2022/isec22_abap/gf_wc.png")



# Upload pentads to GEE ---------------------------------------------------

## Region pentads
ABAP::uploadPentadsToEE(pentads,
                        file.path(ee_get_assethome(),"pentads_wc"),
                        load = FALSE)

## Visit pentads

# Make spatial object and select relevant columns
visits <- gf_abap %>%
    dplyr::left_join(pentads,
                     by = c("Pentad" = "Name")) %>%
    sf::st_sf() %>%
    dplyr::filter(!sf::st_is_empty(.)) %>%     # Remove rows without geometry
    dplyr::mutate(Date = as.character(StartDate)) %>%   # GEE doesn't like dates
    dplyr::select(CardNo, StartDate, Date, Pentad, TotalHours)

# Upload to GEE
visits %>%
    dplyr::select(-c(StartDate, TotalHours)) %>%
    ABAP::uploadPentadsToEE(.,
                            file.path(ee_get_assethome(),"df_visits_wc"),
                            load = FALSE)





# Load pentads from GEE ---------------------------------------------------


# Load pentads from GEE
ee_pentads <- rgee::ee$FeatureCollection(file.path(ee_get_assethome(),"pentads_wc"))

# OR
ee_visits <- rgee::ee$FeatureCollection(file.path(ee_get_assethome(),"df_visits_wc"))



# Annotate with GEE -------------------------------------------------------


## Extract mean NDVI for each pentad
ndvi_mean <- addVarEEcollection(ee_pentads = ee_pentads,
                                collection = "MODIS/006/MOD13A2",
                                dates = c("2010-01-01", "2013-01-01"),
                                temp_reducer = "mean",
                                spt_reducer = "mean",
                                bands = "NDVI")


pentads_terra <- addVarEEcollection(ee_pentads = ee_pentads,                    # Note that we need our remote asset here
                                   collection = "IDAHO_EPSCOR/TERRACLIMATE",   # You can find this in the code snippet
                                   dates = c("2010-01-01", "2011-01-01"),
                                   temp_reducer = "mean",
                                   spt_reducer = "mean",
                                   bands = c("tmmn", "tmmx"))


pentads_water <- addVarEEimage(ee_pentads = ee_pentads,                   # Note that we need our remote asset here
                               image = "JRC/GSW1_3/GlobalSurfaceWater",   # You can find this in the code snippet
                               reducer = "mean",
                               bands = "occurrence",
                               unmask = TRUE)


# Annotate visits with GEE ------------------------------------------------

visits_ndvi <- addVarEEclosestImage(ee_visits,
                                    collection = "MODIS/006/MOD13A2",
                                    reducer = "mean",
                                    maxdiff =  30,
                                    bands = "NDVI")

pentads_water %>%
    ggplot() +
    geom_sf(aes(fill = log(occurrence+0.1))) +
    scale_fill_viridis_c(name = "log water\nrecurrence", direction = -1) +
    theme_void()

ggsave("comms/isec_2022/isec22_abap/water_wc.png")


# Make an unmarked data frame with the data
pentads_water$water2 <- pentads_water$occurrence + 2

gf_um <- abapToUnmarked_single(gf_abap)
summary(gf_um)

gf_um <- addEEtoUnmarked_single(gf_um, pentads_water, c("NDVI"))




# Explore abap pentads ----------------------------------------------------

library(ABAP)
library(sf)
library(dplyr)
library(ggplot2)

url <- "https://api.birdmap.africa/sabap2/v2/pentads?format=geoJSON"

# Extract data
pentads <- sf::read_sf(url)

# Re-classify pentads
sabap2 <- c("South Africa", "eSwatini", "Lesotho", "Namibia", "Zimbabwe", "Zambia", "Botswana", "Mozambique", "Malawi")
pentads <- pentads %>%
    mutate(project = case_when(country %in% sabap2 ~ "SABAP2",
                               country == "Kenya" ~ "Kenya Bird Map",
                               country == "Nigeria" ~ "NiBAP",
                               TRUE ~ "Other"))

plot(pentads["project"], lwd = 0.01)

wc <- rnaturalearth::ne_states(country = 'south africa', returnclass = 'sf') %>%
    filter(name == "Western Cape")

# Plot whole ABAP range
pentads %>%
    # filter(country == "South Africa") %>%
    ggplot() +
    geom_sf(aes(fill = project), lwd = 0) +
    scale_fill_brewer(name = "Project", palette = "Paired",
                      breaks=c('Kenya Bird Map', 'NiBAP', 'SABAP2', 'Other')) +
    geom_sf(data = wc, fill = NA) +
    theme_void()

ggsave(filename = "comms/isec_2022/isec22_abap/abap_map_wc.png")

# Plot whole ABAP range
pentads %>%
    mutate(project = case_when(country %in% sabap2 ~ "SABAP2",
                               TRUE ~ "Other")) %>%
    ggplot() +
    geom_sf(aes(fill = project), lwd = 0) +
    scale_fill_manual(values = c("#af8dc3", "#a6d854"), name = "Project") +
    theme_void()

ggsave(filename = "comms/isec_2022/isec22_abap/sabap_map.png")
