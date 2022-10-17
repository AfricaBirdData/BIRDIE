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


    counts <- utils::read.csv(setSpOutFilePath("abu_model_data", config, sp_code, ".csv"))

    # Create a sequential 'site' variable
    # counts <- counts %>%
    #     dplyr::group_by(LocationCode) %>%
    #     dplyr::mutate(site = dplyr::cur_group_id()) %>%
    #     dplyr::ungroup()

    # Create covariate matrix for summer/winter
    covts_x <- counts %>%
        dplyr::mutate(intcp = 1) %>%
        dplyr::rename_with(~gsub("_mean||_count", "", .x), .cols = dplyr::everything()) %>%
        dplyr::select(intcp, prcp, tmmn, tmmx) %>%
        dplyr::mutate(dplyr::across(.cols = c(prcp, tmmn, tmmx), .fns = ~scale(.x)))

    # Create covariate matrix for mean population change
    covts_u <- counts %>%
        dplyr::filter(season_id == 1) %>%
        dplyr::rename_with(~gsub("_mean||_count", "", .x), .cols = dplyr::everything()) %>%
        dplyr::select(site_id, year_id, pdsi, watext, watrec) %>%
        # dplyr::group_by(site_id) %>%
        # dplyr::mutate(log_diff_pdsi = c(diff(log(pdsi + 1000)), NA)) %>%
        # dplyr::ungroup() %>%
        # dplyr::select(-c(pdsi, diff_pdsi)) %>%
        dplyr::mutate(dplyr::across(.cols = c(pdsi, watext, watrec), .fns = ~scale(.x)))

    # add intercept
    covts_u <- covts_u %>%
        dplyr::mutate(intcp = 1) %>%
        dplyr::select(intcp, dplyr::everything())

    # There is one NA at the end of each time series of covariates because they
    # are lagged variables. Make this zero - it won't affect the results
    covts_u[is.na(covts_u)] <- 0

    # Convert to array of multiple sites
    U <- covts_u %>%
        dplyr::select(-year_id) %>%
        dplyr::nest_by(site_id) %>%
        dplyr::pull(data) %>%
        unlist() %>%
        array(., dim = c(dplyr::n_distinct(counts$year_id),
                         ncol(covts_u) - 2, # this removes year and site columns
                         dplyr::n_distinct(counts$site_id)))


    # Expected (log) summer count for each site would be the mean count over the years
    mean_mu <- counts %>%
        dplyr::filter(season_id == 1) %>%
        dplyr::group_by(site_id) %>%
        dplyr::summarise(mean_count = log(mean(count + 1, na.rm = TRUE))) %>%  # NOTE THE PLUS ONE TO AVOID ZERO COUNTS AND NON-IDENTIFIBILITY OF PARAMETERS
        dplyr::pull(mean_count)


    # Priors for initial states. It shouldn't be much higher than the maximum count at the site
    ini_ub <- counts %>%
        dplyr::group_by(site_id) %>%
        dplyr::summarise(max_count = max(count + 1, na.rm = TRUE)) %>%
        dplyr::pull(max_count) %>%
        log()

    # Prepare data
    data <- list(
        count = counts$count + 1, # NOTE THE PLUS ONE TO AVOID ZERO COUNTS AND NON-IDENTIFIBILITY OF PARAMETERS
        summer = dplyr::case_when(counts$season_id == 1 ~ 1L,
                                  counts$season_id == 2 ~ 0L,
                                  TRUE ~ NA_integer_),
        mean_mu = mean_mu,
        ini_sd = ini_ub - mean_mu,
        nyears = dplyr::n_distinct(counts$year_id),
        nsites = dplyr::n_distinct(counts$site_id),
        year = counts$year_id,
        site = counts$site_id,
        N = nrow(counts),
        X = as.matrix(covts_x),
        U = U,
        K = ncol(covts_x),
        M = ncol(U))

    param = c("phi", "beta", "lambda", "mu.beta", "sig.zeta", "sig.eps", "sig.alpha", "sig.e", "mu_t", "summer")

    # param = c("beta", "mu.beta", "B", "G", "mu_t")

    # Set initial values
    inits <- function(){
        list(ini_s = rnorm(data$nsites, data$mean_mu, data$ini_sd*1.5),
             phi = runif(data$nsites, 0.1, 0.9),
             tau.alpha = rgamma(data$nsites, 3, 2),
             tau.e = rgamma(data$nsites, 3, 2),
             B = matrix(rnorm(data$K * data$nsites, 0, 1), ncol = data$K),
             G = matrix(rnorm(data$M * data$nsites, 0, 1), ncol = data$M),
             tau.eps = rexp(data$nyears-1, 1),
             tau.zeta = rexp(data$nyears-1, 1))
    }

    print(paste("Fitting state-space JAGS model at", Sys.time()))

    # Fit 2-season dynamic trend model
    fit <- jagsUI::jags.basic(data = data,
                              parameters.to.save = param,
                              model.file = "analysis/models/cwac_ssm_lat_season_multi_hier.R",
                              inits = inits,
                              n.chains = 3, n.iter = 10000, n.burnin = 5000,
                              modules = c('glm', 'dic'), parallel = TRUE,
                              n.cores = 3, DIC = TRUE, verbose = TRUE)

    # Save
    if(length(sp_code) > 1){
        sp_code <- "group"
    }

    saveRDS(fit, setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))

}
