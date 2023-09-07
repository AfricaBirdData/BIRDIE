
# Read in site data
site <- readRDS("analysis/data/site_dat_sa_gee_08_19.rds")

# Save as csv
write.csv(sf::st_drop_geometry(site), "analysis/data/site_dat_sa_gee_08_19.csv", row.names = FALSE)


# Read in visit data
visit <- readRDS("analysis/data/visit_dat_sa_gee_08_19.rds")

# Save as csv
write.csv(visit, "analysis/data/visit_dat_sa_gee_08_19.csv", row.names = FALSE)
