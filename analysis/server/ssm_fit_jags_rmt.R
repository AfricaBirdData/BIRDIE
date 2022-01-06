library(BIRDIE)
library(CWAC)
library(dplyr)

rm(list = ls())


# Script settings ---------------------------------------------------------

data_dir <- "/home/birdie/analysis/data/"
mod_file <- "/drv_birdie/Working/git/BIRDIE/analysis/models/cwac_ssm_lat_season.jags"
data_outdir <- "/drv_birdie/birdie_ftp/"
plot_outdir <- "/drv_birdie/birdie_ftp/"


# Prepare count data ------------------------------------------------------

# List all sites in South Africa
# sites <- listCwacSites(.region_type = "country", .region = "South Africa")
#
# # Extract site of interest
# site_spt <- sites %>%
#     filter(LocationName == "Barberspan")
#
# site_id <- site_spt$LocationCode

site_id <- "26352535"

# Years of interest
year_ch <- c("93", "18")

counts <- readRDS(paste0(data_dir, site_id, "_", paste(year_ch, collapse = "_"), "_visit_covts.rds"))

# Identify site
site_id <- counts %>%
    filter(!is.na(LocationCode)) %>%
    pull(LocationCode) %>%
    unique()

if(length(site_id) != 1){
    stop("Either location code is missing or there are multiple sites.")
}

# Get species list
spp <- unique(counts$SppRef)


# Fit models and save results ---------------------------------------------

for(i in seq_along(spp)){

    sp <- spp[i]

    print(paste0("Working on species ", sp, " (", i, " of ", length(spp), ")"))

    # Prepare data to fit an SSM
    ssmcounts <- BIRDIE::prepSsmData(counts, sp, keep = Hmisc::Cs(CountCondition, prcp, tmmn, tmmx, watext, watrec))


    # Prepare covariates ---------------------------------------------------------

    # Create covariate matrix
    covts_x <- ssmcounts %>%
        mutate(intcp = 1) %>%
        dplyr::select(intcp, prcp, tmmn, tmmx, watext, watrec) %>%
        mutate(across(.cols = c(prcp, tmmn, tmmx, watext, watrec), .fns = ~scale(.x)))

    covts_z <- ssmcounts %>%
        mutate(intcp = 1,
               CountCondition = if_else(is.na(CountCondition), 0L, CountCondition)) %>%
        dplyr::select(intcp, CountCondition)


    # JAGS model --------------------------------------------------------------

    # Prepare data (note the addition of 0.1 to avoid infinite values)
    data <- list(obs = log(ssmcounts$count + 0.1),
                 summer = case_when(ssmcounts$Season == "S"~ 1L,
                                    ssmcounts$Season == "W" ~ 0L,
                                    TRUE ~ NA_integer_),
                 nyears = n_distinct(ssmcounts$year),
                 year = ssmcounts %>% group_by(year) %>% mutate(y = cur_group_id()) %>% pull(y),
                 N = nrow(ssmcounts),
                 X = as.matrix(covts_x),
                 K = ncol(covts_x))

    param = c("beta", "lambda", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t")

    # Fit 2-season dynamic trend model
    fit_dyn <- jagsUI::jags(data = data,
                            parameters.to.save = param,
                            model.file = model_file,
                            n.chains = 3, n.iter = 10000, n.burnin = 5000,
                            modules = c('glm','lecuyer', 'dic'), parallel = TRUE,
                            n.cores = 3, DIC = TRUE, verbose = TRUE)

    # Plot
    pers_theme <- ggplot2::theme_bw()
    p <- BIRDIE::plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, linear = TRUE,
                            plot_options = list(pers_theme = pers_theme,
                                                colors = c("#71BD5E", "#B590C7")))

    # grid::grid.newpage()
    # grid::grid.draw(p$plot)

    # Extract plot data
    plotdata <- p$data

    # Prepare state data
    stt_df <- plotdata[[1]]

    stt_dfw <- cbind(stt_df %>%
                         dplyr::filter(season == 1) %>%
                         dplyr::select(-season),
                     stt_df %>%
                         dplyr::filter(season == 2) %>%
                         dplyr::select(-c(year, season)))

    names(stt_dfw) <- c("summer.logest", "summer.logest.ci.lower", "summer.logest.ci.upper",
                        "year", "log.summer.count",
                        "winter.logest", "winter.logest.ci.lower", "winter.logest.ci.upper",
                        "log.winter.count")

    # Prepare trend data
    trd_df <- plotdata[[2]]

    names(trd_df) <- c("slope.est", "slope.ci.lower", "slope.ci.upper",
                       "winter.prop.est", "winter.prop.ci.lower", "winter.prop.ci.upper",
                       "year")

    # Combine
    out_df <- cbind(stt_dfw, dplyr::select(trd_df, -year))

    # Add missing columns
    out_df <- out_df %>%
        dplyr::mutate(species = sp,
                      summer.count = round(exp(log.summer.count)),
                      winter.count = round(exp(log.winter.count)))

    # Order columns
    out_df <- out_df %>%
        dplyr::select(species, year, summer.count, winter.count,
                      log.summer.count, log.winter.count,
                      summer.logest, winter.logest, winter.logest.ci.lower,
                      winter.logest.ci.upper,
                      summer.logest.ci.lower, summer.logest.ci.upper,
                      slope.est, slope.ci.lower, slope.ci.upper,
                      winter.prop.est, winter.prop.ci.lower,
                      winter.prop.ci.upper)

    # Export sample
    datadir <- paste0(data_outdir, sp, "/ssm_dat_", site_id, "_", sp, "new.csv")
    utils::write.csv(out_df, datadir, row.names = FALSE)

    plotdir <- paste0(plot_outdir, sp, "/ssm_plot_", site_id, "_", sp, "new.png")
    ggplot2::ggsave(plotdir, plot = p$plot, width = 6.4, height = 5.2)

}
