#' Fit state-space JAGS model
#'
#' @param sp_code SAFRING code of the species of interest
#' @inheritParams ppl_run_pipe_abu1
#'
#' @return
#' @export
#'
#' @examples
ppl_fit_ssm_model <- function(sp_code, site, year, config, ...){

    # Load counts
    counts <- readRDS(file.path(config$data_dir, paste(site, year, "visit_covts.rds", sep = "_")))

    # Filter years
    counts <- counts %>%
        dplyr::filter(Year %in% config$years)

    # Prepare data to fit an SSM
    ssmcounts <- prepSsmData(counts, sp_code, keep = Hmisc::Cs(CountCondition, prcp, tmmn, tmmx, watext, watrec))

    # Prepare covariates ---------------------------------------------------------

    # Create covariate matrix
    covts_x <- ssmcounts %>%
        dplyr::mutate(intcp = 1) %>%
        dplyr::select(intcp, prcp, tmmn, tmmx, watext, watrec) %>%
        dplyr::mutate(across(.cols = c(prcp, tmmn, tmmx, watext, watrec), .fns = ~scale(.x)))

    covts_z <- ssmcounts %>%
        dplyr::mutate(intcp = 1,
                      CountCondition = if_else(is.na(CountCondition), 0L, CountCondition)) %>%
        dplyr::select(intcp, CountCondition)


    # JAGS model --------------------------------------------------------------

    # Prepare data (note the addition of 0.1 to avoid infinite values)
    data <- list(obs = log(ssmcounts$count + 0.1),
                 summer = dplyr::case_when(ssmcounts$Season == "S"~ 1L,
                                    ssmcounts$Season == "W" ~ 0L,
                                    TRUE ~ NA_integer_),
                 nyears = dplyr::n_distinct(ssmcounts$year),
                 year = ssmcounts %>%
                     dplyr::group_by(year) %>%
                     dplyr::mutate(y = dplyr::cur_group_id()) %>%
                     dplyr::pull(y),
                 N = nrow(ssmcounts),
                 X = as.matrix(covts_x),
                 K = ncol(covts_x))

    param = c("beta", "lambda", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t")

    print(paste("Fitting state-space JAGS model at", Sys.time()))

    # Fit 2-season dynamic trend model
    fit <- jagsUI::jags(data = data,
                        parameters.to.save = param,
                        model.file = config$mod_file,
                        n.chains = 3, n.iter = 10000, n.burnin = 5000,
                        modules = c('glm','lecuyer', 'dic'), parallel = TRUE,
                        n.cores = 3, DIC = TRUE, verbose = TRUE)

    # Save
    if(length(sp_code) > 1){
        sp_code <- "group"
    }

    saveRDS(fit, file.path(config$data_outdir, sp_code, paste0("ssm_fit_", site, "_", config$years_ch, "_", sp_code, ".rds")))

}
