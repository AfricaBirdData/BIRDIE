
library(CWAC)
library(sf)
library(dplyr)
library(ggplot2)

sites <- CWAC:::listCwacSites()

url <- "https://pipeline.birdmap.africa/cwac/location?short=1"

# Extract data
sites <- httr::RETRY("GET", url) %>%
    httr::content(as = "text", encoding = "UTF-8") %>%
    rjson::fromJSON() %>%
    CWAC::jsonToTibble() %>%
    readr::type_convert()

sites <- sites %>%
    select(LocationName, X, Y, Country) %>%
    st_as_sf(coords = c("X", "Y"), crs = 4326, remove = FALSE)

plot(st_geometry(sites))

africa <- rnaturalearth::ne_countries(continent = "Africa", returnclass = 'sf')

plot(st_geometry(africa))
plot(st_geometry(sites), col = "red", pch = 19, add = TRUE)

sites %>%
    filter(Y == max(Y)|
           Y == min(Y))

sites <- sites %>%
    filter(Y < 10,
           Y > -35.2)

# All African CWAC sites
ggplot() +
    geom_sf(data = africa) +
    geom_point(data = sites, aes(x = X, y = Y), col = "red") +
    theme_void()

ggsave(filename = "comms/isec_2022/isec22_abap/cwac_africa.png")

# South African CWAC sites
sa <- rbind(rnaturalearth::ne_states(country = "South Africa", returnclass = 'sf') %>%
                select(admin),
            rnaturalearth::ne_countries(scale = 10, country = "Lesotho", returnclass = 'sf') %>%
                select(admin),
            rnaturalearth::ne_countries(scale = 10, country = "Eswatini", returnclass = 'sf') %>%
                select(admin)) %>%
    st_crop(xmin = 10, xmax = 35, ymin = -40, ymax = -20)

sa %>%
    ggplot() +
    geom_sf() +
    geom_point(data = filter(sites, Country == "South Africa", Y < -20),
               aes(x = X, y = Y), col = "red") +
    theme_void()

ggsave(filename = "comms/isec_2022/isec22_abap/cwac_sa.png")

# Try sites' polygons
url <- "https://pipeline.birdmap.africa/cwac/location?format=geoJSON"

sites_sf <- st_read(url)
sites_sf <- st_make_valid(sites_sf)

plot(st_geometry(sites_sf))


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
    scale_fill_brewer(name = "Project", palette = "Accent",
                      breaks=c('Kenya Bird Map', 'NiBAP', 'SABAP2', 'Other')) +
    geom_sf(data = wc, fill = NA) +
    theme_void()

ggsave(filename = "comms/isec_2022/isec22_abap/abap_map.png")

# Plot Western Cape pentads
pentads %>%
    filter(province == "Western Cape") %>%
    ggplot() +
    geom_sf(fill = NA) +
    theme_voitd


