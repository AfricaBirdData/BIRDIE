library(dplyr)
library(ggplot2)
library(gridExtra)
library(ABAP)
library(CWAC)
library(sf)
library(rnaturalearth)

rm(list = ls())

outdir <- "comms/biodv_informatics"

# Background maps ---------------------------------------------------------

# Download South Africa map with Lesotho and Swaziland
ssa_map <- ne_countries(scale = 10, continent = "Africa", returnclass = "sf")
sa_map <- ne_states("South Africa", returnclass = "sf")

# Crop ssa_map
frame <- raster::extent(c(12, 40, -37, -17))
sf_use_s2(FALSE) # There are some problematic vertices to use spherical geometry, so we switch it off
ssa_map <- st_crop(ssa_map, frame)
sa_map <- st_crop(sa_map, frame)


# ABAP effort -------------------------------------------------------------

# Download South African pentads
sa_pentads <- getRegionPentads("country", "south africa")

# Download ABAP data for South Africa (any species would work)
abapdata <- getAbapData(.spp_code = 6, .region_type = "country", .region = "South Africa")

# Subset to SABAP2 (>2007, it started in mid 2007)
abapdata <- abapdata %>%
    mutate(year = lubridate::year(StartDate)) %>%
    filter(year > 2007, !is.na(year))

# Count number of cards per pentad and plot
abapdata %>%
    count(Pentad, CardNo) %>%
    count(Pentad) %>%
    left_join(sa_pentads %>%
                  dplyr::select(Pentad = pentad)) %>%
    st_sf() %>%
    ggplot() +
    geom_sf(aes(fill = log(n)), lwd = NA) +
    geom_sf(data = ssa_map, fill = NA, size = 0.2) +
    geom_sf(data = sa_map, fill = NA, size = 0.2) +
    scale_fill_viridis_c(option = "B", direction = -1) +
    coord_sf(xlim = c(16, 34), ylim = c(-35, -22)) +
    theme_bw() +
    theme(text = element_text(size = 16))

ggsave(filename = file.path(outdir, "figures_paper", "sabap_effort.png"))


# CWAC effort -------------------------------------------------------------

# Download South African sites
sites <- listCwacSites("country", "south africa")

# Create variables for first and last year monitored
sites <- sites %>%
    mutate(lastyear = lubridate::year(LastSurvey),
           firstyear = lubridate::year(FirstSurvey))

# Remove NA from first and last year, and years prior to 1990
sites <- sites %>%
    filter(!is.na(firstyear),
           !is.na(lastyear),
           firstyear > 1990)

# Filter South African sites and remove a site in the ocean
sites <- sites %>%
    filter(Y < -10, Y > -35)

# Calculate duration and filter out variables
sites <- sites %>%
    mutate(duration = lastyear - firstyear) %>%
    arrange(firstyear, duration) %>%
    select(LocationCode, X, Y, firstyear, lastyear, duration) %>%
    distinct() %>%
    mutate(LocationCode = factor(LocationCode, levels = .$LocationCode))

# Plot timeline
sites %>%
    ggplot() +
    geom_segment(aes(x = firstyear, xend = lastyear,
                     y = LocationCode, yend = LocationCode, col = duration), size = 0.5) +
    scale_color_viridis_c(name = "dur (years)") +
    xlab("Year") + ylab("Site") +
    ggtitle("Timeline of site counts (each fine horizontal line is a site)") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_line())

ggsave(filename = file.path(outdir, "figures_paper", "cwac_timeline.png"))

# Create a dataframe that capture the number CWAC sites active, added and removed
# each year
sites_year <- data.frame(active = integer(length = 31),
                         added = integer(length = 31),
                         removed = integer(length = 31),
                         year = integer(length = 31))

for(i in seq_len(nrow(sites_year))){
    yy <- 1990 + i

    sites_year$year[i] <- yy

    sites_year$active[i] <- sites %>%
        mutate(monitored = ifelse(firstyear <= yy & lastyear >= yy, 1, 0)) %>%
        pull(monitored) %>%
        sum()

    sites_year$added[i] <- sites %>%
        filter(firstyear == yy) %>%
        nrow()

    sites_year$removed[i] <- sites %>%
        filter(lastyear == yy) %>%
        nrow()

}

col_sel <- RColorBrewer::brewer.pal(3, "Dark2")

refact <- 3 # for secondary axis
sites_per_year <- sites_year %>%
    tidyr::pivot_longer(-year, names_to = "state", values_to = "value") %>%
    mutate(value = ifelse(state == "active", value, value * refact)) %>%
    ggplot(aes(x = year, y = value, col = state, group = state)) +
    geom_point() +
    geom_line() +
    scale_color_brewer(name = "", palette = "Dark2", labels = c("Active", "First", "Last")) +
    scale_y_continuous("Sites active",
                       sec.axis = sec_axis(~ . /refact, name = "Sites first/last")) +
    xlab("Year") +
    theme_bw() +
    theme(text = element_text(size = 12))

ggsave(sites_per_year, filename = file.path(outdir, "figures_paper", "number_cwac_sites.png"))


# make separate plots for active and added/removed
active <- sites_year %>%
    tidyr::pivot_longer(-year, names_to = "state", values_to = "value") %>%
    filter(state == "active") %>%
    ggplot(aes(x = year, y = value), col = col_sel[1]) +
    geom_point() +
    geom_line() +
    xlab("Year") + ylab("Number of sites") +
    theme_bw() +
    theme(text = element_text(size = 16))

add_rem <- sites_year %>%
    tidyr::pivot_longer(-year, names_to = "state", values_to = "value") %>%
    filter(state != "active") %>%
    ggplot(aes(x = year, y = value, col = state, group = state)) +
    geom_point() +
    geom_line() +
    scale_color_manual(name = "", values = col_sel[-1], labels = c("First", "Last")) +
    xlab("Year") + ylab("Number of sites") +
    theme_bw() +
    theme(text = element_text(size = 16))


# Plot sites spatially
sites_spt <- sites %>%
    st_as_sf(coords = c("X", "Y"), dim = "XY", remove = FALSE, crs = 4326)

distrib <- sites_spt %>%
    ggplot() +
    geom_sf(aes(fill = duration), shape = 21, stroke = 0.2) +
    geom_sf(data = ssa_map, fill = NA, size = 0.2) +
    geom_sf(data = sa_map, fill = NA, size = 0.2) +
    scale_fill_viridis_c(name = "Years", option = "B", direction = -1) +
    coord_sf(xlim = c(16, 34), ylim = c(-35, -22)) +
    theme_bw() +
    theme(text = element_text(size = 12))

ggsave(filename = file.path(outdir, "figures_paper", "cwac_map.png"))

# Combine into a single figure
ly_mat <- rbind(c(NA, 1, NA),
                c(2, 2, 2),
                c(2, 2, 2))

effort_plot <- arrangeGrob(sites_per_year, distrib,
                           layout_matrix = ly_mat,
                           widths = c(0.11, 0.8, 0.04),
                           heights = c(1/3, 1/3, 1/3))


ggsave(effort_plot, filename = file.path(outdir, "figures_paper", "cwac_effort.png"))

