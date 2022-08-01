library(BIRDIE)

rm(list = ls())

config <- BIRDIE::configPreambJAGS(2017, server = TRUE)

rhat_df <- data.frame()

for(i in 1:length(config$species)){

    sp_code <- config$species[i]

    print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

    # Skip species if less than 5 suitable sites were detected during fitting model fitting
    error_file <- setSpOutFilePath("Less_5_sites", config, sp_code, ".txt")

    if(file.exists(error_file)){
        message(paste0("Less_5_sites_", sp_code, "_", config$years_ch, ".txt"))
        next
    }

    # Else proceed with diagnostics
    fit <- readRDS(setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))

    fit_stats <- BIRDIE:::processJAGSoutput(fit, DIC = FALSE, params.omit = NULL)

    rhats <- lapply(fit_stats$Rhat, function(x) sum(abs(x - 1) > 0.1)) %>%
        as.data.frame()

    rhats$sp <- sp_code
    rhats$nobs <- length(fit_stats$mean$mu_t)
    rhats$years <- config$years_ch

    rhat_df <- rbind(rhat_df, rhats)

}

write.csv(rhat_df,
          file.path(config$out_dir, paste0("abu_diag_", config$years_ch, ".csv")),
          row.names = FALSE)
