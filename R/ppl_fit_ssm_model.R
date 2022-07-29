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
    counts <- counts %>%
        dplyr::group_by(LocationCode) %>%
        dplyr::mutate(site = dplyr::cur_group_id()) %>%
        dplyr::ungroup()

    # Create covariate matrix for summer/winter
    covts_x <- counts %>%
        dplyr::mutate(intcp = 1) %>%
        dplyr::rename_with(~gsub("_mean||_count", "", .x), .cols = dplyr::everything()) %>%
        dplyr::select(intcp, prcp, tmmn, tmmx) %>%
        dplyr::mutate(dplyr::across(.cols = c(prcp, tmmn, tmmx), .fns = ~scale(.x)))

    # Create covariate matrix for mean population change
    covts_u <- counts %>%
        dplyr::filter(Season == "S") %>%
        dplyr::rename_with(~gsub("_mean||_count", "", .x), .cols = dplyr::everything()) %>%
        dplyr::select(year, site, pdsi, watext, watrec) %>%
        dplyr::group_by(year, site) %>%
        dplyr::summarize(dplyr::across(dplyr::everything(), list(mean))) %>%
        dplyr::rename_with(~gsub("_1", "", .x), .cols = dplyr::everything()) %>%
        dplyr::ungroup() %>%
        # dplyr::mutate(intcp = 1) %>%
        dplyr::group_by(site) %>%
        dplyr::arrange(site, year) %>%
        dplyr::mutate(dplyr::across(.cols = c(pdsi, watext, watrec), .fns = ~dplyr::lead(.x) - .x)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(dplyr::across(.cols = c(pdsi, watext, watrec), .fns = ~scale(.x)))

    # Convert to array of multiple sites
    U <- covts_u %>%
        dplyr::select(-year) %>%
        dplyr::nest_by(site) %>%
        dplyr::pull(data) %>%
        unlist() %>%
        array(., dim = c(dplyr::n_distinct(counts$year),
                         ncol(covts_u) - 2, # this removes year and site columns
                         dplyr::n_distinct(counts$site)))

    # Prepare data
    data <- list(
        count = counts$count + 1, # NOTE THE PLUS ONE TO AVOID ZERO COUNTS AND NON-IDENTIFIBILITY OF PARAMETERS
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
        U = U,
        K = ncol(covts_x),
        M = ncol(U))

    param = c("phi", "beta", "lambda", "mu.beta", "sig.zeta", "sig.eps", "sig.alpha", "sig.e", "mu_t")

    # param = c("beta", "mu.beta", "B", "G", "mu_t")

    # Set initial values
    inits <- function(){
        list(phi = runif(data$nsites, 0, 1),
             tau.alpha = rgamma(data$nsites, 2, 2),
             tau.e = rgamma(data$nsites, 2, 2),
             B = matrix(rnorm(data$K * data$nsites, 0, 2), ncol = data$K),
             G = matrix(rnorm(data$M * data$nsites, 0, 2), ncol = data$M),
             tau.eps = rexp(data$nyears-1, 0.5),
             tau.zeta = rexp(data$nyears-1, 0.5))
    }

    print(paste("Fitting state-space JAGS model at", Sys.time()))

    # Fit 2-season dynamic trend model
    fit <- jagsUI::jags.basic(data = data,
                              parameters.to.save = param,
                              model.file = "analysis/models/cwac_ssm_lat_season_multi_hier.R",
                              inits = inits,
                              n.chains = 3, n.iter = 15000, n.burnin = 10000,
                              modules = c('glm', 'dic'), parallel = TRUE,
                              n.cores = 3, DIC = TRUE, verbose = TRUE)

    # Save
    if(length(sp_code) > 1){
        sp_code <- "group"
    }

    saveRDS(fit, setSpOutFilePath("ssm_fit", config, sp_code, ".rds"))

}
