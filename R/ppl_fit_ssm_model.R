#' Fit state-space JAGS model
#'
#' @param sp_code SAFRING code of the species of interest
#' @inheritParams ppl_run_pipe_abu1
#'
#' @return
#' @export
#'
#' @examples
ppl_fit_ssm_model <- function(sp_code, config, ...){


    counts <- read.csv(setSpOutFilePath("abu_model_data", config, sp_code, ".csv"))

    # Create a sequential 'site' variable
    counts <- counts %>%
        dplyr::group_by(LocationCode) %>%
        dplyr::mutate(site = dplyr::cur_group_id()) %>%
        dplyr::ungroup()

    # Create covariate matrix
    covts_x <- counts %>%
        dplyr::mutate(intcp = 1) %>%
        dplyr::rename_with(~gsub("_mean||_count", "", .x), .cols = dplyr::everything()) %>%
        dplyr::select(intcp, prcp, tmmn, tmmx, pdsi, watext, watrec) %>%
        dplyr::mutate(dplyr::across(.cols = c(prcp, tmmn, tmmx, pdsi, watext, watrec), .fns = ~scale(.x)))

    # Prepare data
    data <- list(
        count = counts$count,
        summer = dplyr::case_when(counts$Season == "S"~ 1L,
                                  counts$Season == "W" ~ 0L,
                                  TRUE ~ NA_integer_),
        nyears = dplyr::n_distinct(counts$year),
        nsites = dplyr::n_distinct(counts$site),
        year = counts %>%
            dplyr::group_by(year) %>%
            dplyr::mutate(y = dplyr::cur_group_id()) %>%
            dplyr::pull(y),
        site = counts$site,
        N = nrow(counts),
        X = as.matrix(covts_x),
        K = ncol(covts_x))

    param = c("beta", "lambda", "sig.zeta", "sig.eps", "sig.alpha", "sig.e", "mu_t")

    print(paste("Fitting state-space JAGS model at", Sys.time()))

    # Fit 2-season dynamic trend model
    fit <- jagsUI::jags(data = data,
                        parameters.to.save = param,
                        model.file = "analysis/models/cwac_ssm_lat_season_multi_hier.R",
                        n.chains = 3, n.iter = 15000, n.burnin = 10000,
                        modules = c('glm', 'dic'), parallel = TRUE,
                        n.cores = 3, DIC = TRUE, verbose = TRUE)

    # Save
    if(length(sp_code) > 1){
        sp_code <- "group"
    }

    saveRDS(fit, setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))

}
