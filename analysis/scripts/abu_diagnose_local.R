library(BIRDIE)

rm(list = ls())

config <- BIRDIE::configPreambJAGS(2017, server = FALSE)

rhat_df <- data.frame()

for(i in 1:length(config$species)){

    sp_code <- config$species[i]

    print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

    # Skip species if less than 5 suitable sites were detected during fitting model fitting
    error_file <- file.path(config$out_dir, sp_code, paste0("Less_5_sites_", sp_code, "_", config$years_ch, ".txt"))

    if(file.exists(error_file)){
        message(paste0("Less_5_sites_", sp_code, "_", config$years_ch, ".txt"))
        next
    }

    # Else proceed with diagnostics
    fit <- readRDS(file.path(config$out_dir, sp_code, paste0("ssm_fit_", config$years_ch, "_", sp_code, ".rds")))

    names(fit$Rhat)

    rhats <- lapply(fit$Rhat, function(x) sum(x > 1.09)) %>%
        as.data.frame()

    rhats$sp <- sp_code
    rhats$nobs <- length(fit$mean$mu_t)
    rhats$years <- config$years_ch

    rhat_df <- rbind(rhat_df, rhats)

}

write.csv(rhat_df,
          file.path(config$out_dir, paste0("abu_diag_", config$years_ch, "_", sp_code, ".csv")),
          row.names = FALSE)
