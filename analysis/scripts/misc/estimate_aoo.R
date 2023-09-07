library(dplyr)
library(BIRDIE)
library(occuR)

rm(list = ls())

data_dir <- "analysis/data/"
fit_dir <- "analysis/out_nosync/"

sp_sel <- 235
years <- 2008

sp_detect <- ABAP::getAbapData(.spp_code = sp_sel,
                               .region_type = "country",
                               .region = "South Africa",
                               .years = years)


# Load occupancy data -----------------------------------------------------

# Load visit data, subset years and add detections
visitdata <- readRDS(paste0(data_dir, "visit_dat_sa_gee_08_19.rds")) %>%
    filter(year %in% years) %>%
    left_join(sp_detect %>%
                  dplyr::select(CardNo, obs = Spp) %>%
                  mutate(obs = if_else(obs == "-", 0, 1)),
              by = "CardNo")

# Load site data
sitedata <- readRDS(paste0(data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    dplyr::select(Name, lon, lat, watocc_ever, dist_coast, ends_with(match = as.character(years)))

# I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE
sitedata <- sitedata %>%
    tidyr::drop_na()


# Format to occuR ---------------------------------------------------------

occuRdata <- prepDataOccuR(site_data = sitedata %>%
                               sf::st_drop_geometry() %>%
                               gatherYearFromVars(vars = names(.)[-c(1:5)], sep = "_") %>% # check that 1:5 are the variables that don't change over time
                               mutate(tdiff = tmmx - tmmn),
                           visit_data = visitdata,
                           scaling = list(visit = NULL,
                                          site = c("dist_coast", "prcp", "tdiff", "ndvi", "watext", "watrec")))

# Remove data from missing pentads? (THIS SHOULDN'T HAPPEN WHEN PENTADS ARE DOWNLOADED FROM THE API)
occuRdata$visit <- occuRdata$visit %>%
    dplyr::filter(!is.na(site))

occuRdata$site <- occuRdata$site %>%
    dplyr::filter(site %in% unique(occuRdata$visit$site))

if(years < 2010) year_ch <- "08_10"

# Fit
fit <- readRDS(paste0(fit_dir, sp_sel, "/occur_fit_", year_ch, "_", sp_sel, ".rds"))

# Predict
pred_data <- prepPredictDataOccuR(occuRdata, sf::st_drop_geometry(sitedata),
                                  years = years, scaling = TRUE)

pred <- predict(fit, occuRdata$visit, pred_data, nboot = 1000)


# Estimate realized occupancy ---------------------------------------------

real_occu <- summarizePredOccuR(pred_p = as.vector(pred$p),
                                pred_psi = as.vector(pred$psi),
                                pred_data = pred_data,
                                visit_data = occuRdata$visit)

# Plot
library(ggplot2)
sitedata %>%
    left_join(real_occu, by = "Name") %>%
    ggplot() +
    geom_sf(aes(fill = real_occu), lwd = 0.01) +
    scale_fill_viridis_c()

occu <- matrix(nrow = 1000, ncol = nrow(sitedata))

st <- Sys.time()
for(i in 1:1000){

    pred_occu <- summarizePredOccuR(pred_p = pred$pboot[i,],
                                    pred_psi = pred$psiboot[i,],
                                    pred_data = pred_data,
                                    visit_data = occuRdata$visit)

    occu[i,] <- pred_occu$real_occu

}

en <- Sys.time()
en - st

aoo <- rowSums(occu)

hist(aoo, freq = F)
summary(aoo);sd(aoo)
curve(dnorm(x, 437.7, 23.85), add = T)
