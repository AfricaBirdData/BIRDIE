# In this script we annotate CWAC data with environmental data from Google
# Earth Engine. We use the site Barberspan as an example, but the script
# should be adapted to run all CWAC sites.

library(CWAC)
library(dplyr)
library(sf)
library(ABAP)
library(rgee)

ee_check()
ee_Initialize(drive = TRUE)

rm(list = ls())

data_dir <- "analysis/data/"
data_outdir <- "analysis/data/"


# Prepare site data -------------------------------------------------------

# List all sites in South Africa
sites <- listCwacSites(.region_type = "country", .region = "South Africa")

# Extract site of interest
site_spt <- sites %>%
  filter(LocationName == "Barberspan") %>%
  select(LocationCode, LocationName, X, Y) %>%
  st_as_sf(coords = c("X", "Y"), dim = "XY", crs = st_crs(4326))

site_id <- site_spt$LocationCode

# Read in catchment data
catchmt <- read_sf(paste0(data_dir, "catchmt_4.shp"))

# Extract catchment
site_ctm <- catchmt[unlist(st_intersects(site_spt, catchmt)),]


# Get counts for the site -------------------------------------------------

# Get counts between 93 and 20
year_rg <- c(1993, 2020)
counts <- getCwacSiteCounts(site_id)
counts <- counts %>%
  filter(Year >= year_rg[1], Year <= year_rg[2])

# Add missing surveys
counts <- addMissingCwacCounts(counts, years = year_rg[1]:year_rg[2])

# Give missing surveys a date based on the dates from other surveys
month_summer <- counts %>%
  mutate(month = lubridate::month(StartDate)) %>%
  filter(Season == "S", !is.na(month)) %>%
  count(month) %>%
  filter(n == max(n)) %>%
  pull(month)

month_winter <- counts %>%
  mutate(month = lubridate::month(StartDate)) %>%
  filter(Season == "W", !is.na(month)) %>%
  count(month) %>%
  filter(n == max(n)) %>%
  pull(month)

counts <- counts %>%
  mutate(StartDate = case_when(is.na(StartDate) & Season == "S" ~ as.Date(paste(Year, month_summer, "01", sep = "-")),
                               is.na(StartDate) & Season == "W" ~ as.Date(paste(Year, month_winter, "01", sep = "-")),
                               TRUE ~ StartDate) ) %>%
  arrange(StartDate, Season)

# Extract visit dates, add geometry and id
visits <- counts %>%
  distinct(Card, StartDate) %>%
  mutate(id = row_number(),
         geometry = site_ctm$geometry) %>%
  st_sf()

# Find year range
year_rg <- range(lubridate::year(visits$StartDate))
year_ch <- substr(year_rg, nchar(year_rg) - 1, nchar(year_rg))


# Create output name
outfile <- paste0(data_outdir, site_id, "_", paste(year_ch, collapse = "_"), "_visit_covts.rds")

# Annotate with TerraClimate ----------------------------------------------

# Set a name for the asset
eeid <- sprintf("%s/%s", ee_get_assethome(), 'cwac_visits')

# Upload to EE (if not done already)
visits %>%
  rename(Date = StartDate) %>%
  mutate(Date = lubridate::floor_date(Date, "month")) %>%
  mutate(Date = as.character(Date)) %>%
  uploadPentadsToEE(asset_id = eeid,
                    load = FALSE)

# Load the remote data asset
ee_visit <- ee$FeatureCollection(eeid)

# Define bands
bands <- c("pr", "tmmn", "tmmx")

# Define function
f <- function(band){

  # Annotate with GEE TerraClimate
  visit_env <- addVarEEclosestImage(ee_pentads = ee_visit,
                                    collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                    reducer = "mean",
                                    maxdiff = 15,
                                    bands = band)

  # Fix names and variables
  visit_env <- visit_env %>%
    rename_with(~gsub("val", band, .x), .cols = starts_with("val")) %>%
    dplyr::select(id, all_of(band)) %>%
    st_drop_geometry()

  return(visit_env)

}

# Annotate
visit_vars <- purrr::map(bands, ~f(.x))

# bind
visit_vars <- visits %>%
  st_drop_geometry() %>%
  left_join(bind_cols(visit_vars) %>%
              rename(id = id...1) %>%
              dplyr::select(id, all_of(bands)),
            by = "id")

# Fix covariates
visit_vars <- visit_vars %>%
  dplyr::select(id, everything()) %>%
  rename_with(.fn =  ~gsub("pr", "prcp", .x),
              .cols = starts_with("pr")) %>%
  rename_with(.fn =  ~gsub("tmmn", "tmmn", .x),
              .cols = starts_with("tmmn")) %>%
  rename_with(.fn =  ~gsub("tmmx", "tmmx", .x),
              .cols = starts_with("tmmx")) %>%
  mutate(across(.cols = starts_with("tmmn"),
                .fns = ~.x/10),
         across(.cols = starts_with("tmmx"),
                .fns = ~.x/10))

saveRDS(visit_vars, outfile)


# Annotate with water occurrence ------------------------------------------

# Set a name for the asset
eeid <- sprintf("%s/%s", ee_get_assethome(), 'cwac_visits')

# Upload to EE (if not done already)
visits %>%
  rename(Date = StartDate) %>%
  mutate(Date = lubridate::floor_date(Date, "year")) %>%
  mutate(Date = as.character(Date)) %>%
  uploadPentadsToEE(asset_id = eeid,
                    load = FALSE)

# Load the remote data asset
ee_visit <- ee$FeatureCollection(eeid)

visit_water <- vector("list", length = 2)

# Number of pixels with water each year
visit_water[[1]] <- addVarEEclosestImage(ee_pentads = ee_visit,
                                         collection = "JRC/GSW1_3/YearlyHistory",
                                         reducer = "count",
                                         maxdiff = 200,
                                         bands = "waterClass")

# Fix names and variables
visit_water[[1]] <- visit_water[[1]] %>%
  rename_with(~gsub("val", "watext", .x), .cols = starts_with("val")) %>%
  dplyr::select(id, all_of("watext")) %>%
  st_drop_geometry()

# Recurrence of pixels with water each year
visit_water[[2]] <- addVarEEclosestImage(ee_pentads = ee_visit,
                                         collection = "JRC/GSW1_3/YearlyHistory",
                                         reducer = "mean",
                                         maxdiff = 200,
                                         bands = "waterClass")

# Fix names and variables
visit_water[[2]] <- visit_water[[2]] %>%
  rename_with(~gsub("val", "watrec", .x), .cols = starts_with("val")) %>%
  dplyr::select(id, all_of("watrec")) %>%
  st_drop_geometry()

# bind
visit_vars <- visit_water[[1]] %>%
  left_join(bind_cols(visit_water[[2]]),
            by = "id")

# join with previous variables and save
readRDS(outfile) %>%
  left_join(visit_vars, by = "id") %>%
  saveRDS(outfile)



# Save data with covariates -----------------------------------------------

visit_vars <- readRDS(outfile)

counts %>%
  left_join(dplyr::select(visit_vars, -id), by = c("Card", "StartDate")) %>%
  saveRDS(outfile)
