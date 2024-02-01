#' Add detection covariate to spOccupancy data list
#'
#' @param spOcc_data An spOccupancy data list
#' @param covt_data A data frame with the ABAP visits used to create spOcc_data.
#' The data frame must contain the columns 'Pentad', 'StartDate' and the covariate
#' that we want to add to spOcc_data. If the data is multi-season, then we also
#' need the variable that identifies the season. Only one covariate at a time
#' is allowed at the moment.
#' @param seasons The name of the variable used to identify the seasons in
#' multi-season data sets. Defaults to NULL.
#'
#' @return An spOccupancy data list with the additional detection covariate.
#' @export
#'
#' @examples
addSpOccDetCovt <- function(spOcc_data, covt_data, seasons = NULL){

    # Sort data frame for consistency with other functions
    covt_data <- covt_data %>%
        dplyr::arrange(Pentad, StartDate)

    # Extract unique pentads
    pentad_id <- dimnames(spOcc_data$y)[[1]]
    n_sites <- length(pentad_id)

    # Extract maximum number of visits in a single season
    if(is.null(seasons)){

        max_visits <- dim(spOcc_data$y)[2]
        n_seasons <- 1
        season_vec <- 1

        covt_data <- covt_data %>%
            dplyr::mutate(season = as.character(1))

    } else {

        if(!seasons %in% names(covt_data)){
            stop(paste("The field", seasons, "specified as season, is not present in the data"))
        }

        # Extract seasons and calculate max number of visits per season
        # Note that we will need to match to dimnames of spOcc_data with are characters
        covt_data <- covt_data %>%
            dplyr::mutate(season = !!as.name(seasons)) %>%
            dplyr::mutate(season = as.character(season))

        max_visits <- dim(spOcc_data$y)[3]
        n_seasons <- dim(spOcc_data$y)[2]
        season_vec <- dimnames(spOcc_data$y)[[2]]

        # Sometimes there might not be names. We make new default ones
        if(is.null(season_vec)){
            season_vec <- seq_len(n_seasons)
        }

    }

    # Aux padding vector
    vpad <- rep(NA, max_visits)

    # Aux pentad column to join yearly results
    pentad_col <- data.frame(Pentad = pentad_id)

    # Extract variable name
    var_names <- names(covt_data)[which(!names(covt_data) %in% c("Pentad", "StartDate", "season", seasons))]

    for(v in seq_along(var_names)){

        var_name <- var_names[v]

        ## Create empty array
        covt_out <-  array(NA, dim = c(n_sites, n_seasons, max_visits),
                           dimnames = list(pentad_id, season_vec, seq_len(max_visits)))

        for(k in seq_along(season_vec)){

            ## Create dataframe to format
            season_data <- covt_data %>%
                dplyr::filter(season == season_vec[k]) %>%
                # dplyr::select(Pentad, Spp, TotalHours, StartDate) %>%
                # dplyr::mutate(Spp = ifelse(Spp == "-", 0L, 1L),
                #               julian_day = lubridate::yday(StartDate)) %>%
                dplyr::right_join(pentad_col, by = "Pentad") %>%      # Join in Pentads missing for the year
                dplyr::nest_by(Pentad) %>%
                dplyr::arrange(Pentad) %>%
                dplyr::mutate(varpad = list(head(c(data[, var_name, drop = TRUE], vpad), max_visits)))

            ## Extract total hours
            covt_out[,k,] <- do.call("rbind", season_data$varpad)

        }

        # If there are no seasons then return a matrix
        if(length(season_vec) == 1){

            covt_out <- as.matrix(covt_out[, 1, ])

        }


        # Make data list
        if(is.null(spOcc_data$det.covs)) {
            spOcc_data$det.covs <- covt_out
        } else {
            spOcc_data$det.covs <- c(spOcc_data$det.covs, list(covt_out))
            names(spOcc_data$det.covs)[length(spOcc_data$det.covs)] <- var_name
        }
    }
    return(spOcc_data)

}



defineSpOccuPriors <- function(prev_fit, config){

    # The previous model might not be exactly the same as the new model so the
    # priors need to match the model covariates

    # We first set generic priors and then substitute those that we have
    # previous estimates for

    # Extract fixed effects from config
    fx_occ <- config$occ_mod[!grepl("\\|", config$occ_mod)]
    fx_det <- config$det_mod[!grepl("\\|", config$det_mod)]

    # add intercepts?
    fx_occ <- c("(Intercept)", fx_occ)
    fx_det <- c("(Intercept)", fx_det)

    nocc <- length(fx_occ)
    ndet <- length(fx_det)

    prior_mean_alpha <- rep(0, ndet)
    prior_mean_beta <- rep(0, nocc)
    prior_var_alpha <- rep(2.5, ndet)
    prior_var_beta <- rep(2.5, nocc)

    names(prior_mean_alpha) <- names(prior_var_alpha) <- fx_det
    names(prior_mean_beta) <- names(prior_var_beta) <- fx_occ

    # Extract new priors from fit
    new_prior_mean_beta <- apply(prev_fit$beta.samples, 2, mean)
    new_prior_sd_beta <- apply(prev_fit$beta.samples, 2, sd) + 1
    new_prior_sd_beta <- pmin(new_prior_sd_beta, sqrt(2.5))  # Set maximum sd to the default
    # names(new_prior_mean_beta) <- names(new_prior_sd_beta) <- dimnames(prev_fit$beta.samples)[[2]]

    new_prior_mean_alpha <- apply(prev_fit$alpha.samples, 2, mean)
    new_prior_sd_alpha <- apply(prev_fit$alpha.samples, 2, sd) + 1
    new_prior_sd_alpha <- pmin(new_prior_sd_alpha, sqrt(2.5))   # Set maximum sd to the default
    # names(new_prior_mean_alpha) <- names(new_prior_sd_alpha) <- dimnames(prev_fit$alpha.samples)[[2]]

    # We want the intercepts to be freely estimated for each year (these are the default priors)
    new_prior_mean_beta["(Intercept)"] <- new_prior_mean_alpha["(Intercept)"] <- 0
    new_prior_sd_beta["(Intercept)"] <- new_prior_sd_alpha["(Intercept)"] <- sqrt(2.5)

    # Substitute where necessary
    prior_mean_alpha[names(prior_mean_alpha) %in% names(new_prior_mean_alpha)] <- new_prior_mean_alpha[match(names(prior_mean_alpha), names(new_prior_mean_alpha))]
    prior_mean_beta[names(prior_mean_beta) %in% names(new_prior_mean_beta)] <- new_prior_mean_beta[match(names(prior_mean_beta), names(new_prior_mean_beta))]
    prior_var_alpha[names(prior_var_alpha) %in% names(new_prior_sd_alpha)] <- new_prior_sd_alpha[match(names(prior_var_alpha), names(new_prior_sd_alpha))]^2
    prior_var_beta[names(prior_var_beta) %in% names(new_prior_sd_beta)] <- new_prior_sd_beta[match(names(prior_var_beta), names(new_prior_sd_beta))]^2

    priors <- list(alpha.normal = list(mean = prior_mean_alpha,
                                       var = prior_var_alpha),
                   beta.normal = list(mean = prior_mean_beta,
                                      var = prior_var_beta))


    # Add random effects if necessary
    re_occ <- config$occ_mod[grepl("\\|", config$occ_mod)]
    re_det <- config$det_mod[grepl("\\|", config$det_mod)]

    re_occ <- gsub("\\(|\\||\\)|1", "", re_occ)
    re_det <- gsub("\\(|\\||\\)|1", "", re_det)

    if(length(re_occ) != 0){
        shape.psi.ig <- rep(1, length(re_occ))
        scale.psi.ig <- rep(1.5, length(re_occ))
        names(shape.p.ig) <- names(scale.psi.ig) <- re_occ
    } else {
        shape.psi.ig <- NULL
        scale.psi.ig <- NULL
    }

    if(length(re_det) != 0){
        shape.p.ig <- rep(1, length(re_det))
        scale.p.ig <- rep(1.5, length(re_det))
        names(shape.p.ig) <- names(scale.p.ig) <- re_det
    } else {
        shape.p.ig <- NULL
        scale.p.ig <- NULL
    }


    # Extract random effects from previous fit
    if(!is.null(prev_fit$sigma.sq.p.samples)){
        new_prior_mean_sigma.sq.p <- apply(prev_fit$sigma.sq.p.samples, 2, mean)
        new_prior_sd_sigma.sq.p <- apply(prev_fit$sigma.sq.p.samples, 2, sd) + 1
        new_prior_sd_sigma.sq.p <- pmin(new_prior_sd_sigma.sq.p, sqrt(2.5))   # Set maximum sd to the default
        # names(new_prior_mean_sigma.sq.p) <- names(new_prior_sd_sigma.sq.p) <- dimnames(prev_fit$sigma.sq.p.samples)[[2]]

        new_scale.p.ig <- new_prior_sd_sigma.sq.p^2 / new_prior_mean_sigma.sq.p
        new_shape.p.ig <- new_prior_mean_sigma.sq.p / new_scale.p.ig

        # Substitute where necessary
        keep <- match(names(scale.p.ig), names(new_scale.p.ig))
        keep <- keep[!is.na(keep)]

        scale.p.ig[which(names(scale.p.ig) %in% names(new_scale.p.ig))] <- new_scale.p.ig[keep]
        shape.p.ig[which(names(shape.p.ig) %in% names(new_shape.p.ig))] <- new_shape.p.ig[keep]

    }

    if(!is.null(prev_fit$sigma.sq.psi.samples)){
        new_prior_mean_sigma.sq.psi <- apply(prev_fit$sigma.sq.psi.samples, 2, mean)
        new_prior_sd_sigma.sq.psi <- apply(prev_fit$sigma.sq.psi.samples, 2, sd) + 1
        new_prior_sd_sigma.sq.psi <- pmin(new_prior_sd_sigma.sq.psi, sqrt(2.5))   # Set maximum sd to the default
        # names(new_prior_mean_sigma.sq.psi) <- names(new_prior_sd_sigma.sq.psi) <- dimnames(prev_fit$sigma.sq.psi.samples)[[2]]

        new_scale.psi.ig <- new_prior_sd_sigma.sq.psi^2 / new_prior_mean_sigma.sq.psi
        new_shape.psi.ig <- new_prior_mean_sigma.sq.psi / new_scale.psi.ig

        # Substitute where necessary
        keep <- match(names(scale.psi.ig), names(new_scale.psi.ig))
        keep <- keep[!is.na(keep)]

        scale.psi.ig[which(names(scale.psi.ig) %in% names(new_scale.psi.ig))] <- new_scale.psi.ig[keep]
        shape.psi.ig[which(names(shape.psi.ig) %in% names(new_shape.psi.ig))] <- new_shape.psi.ig[keep]

    }


    # Add random effects
    priors <- c(priors,
                list(sigma.sq.psi.ig = list(shape = shape.psi.ig,
                                            scale = scale.psi.ig)))

    priors <- c(priors,
                list(sigma.sq.p.ig = list(shape = shape.p.ig,
                                          scale = scale.p.ig)))

    return(priors)

}



#' Conduct goodness-of-fit test for spOccupancy models
#' @description This function comes largely from the ppcOcc.R from the
#' \href{https://github.com/doserjef/spOccupancy}{spOccupancy} package
#' @param object an spOccupancy fit
#' @param post_sims A list of posterior simulations obtained with simDetSpOccu().
#' @param fit_stat Goodness of fit statistic to compute. Currently only Chi-squared
#' is supported.
#' @param group Whether the GOF statistic should be computed for sites (1) or
#' for visits (2). Currently, only grouping by site is supported
#'
#' @return A list GOF statistics computed for observed data and for simulated
#' data. One statistic for is computed for each MCMC iteration.
#' The function will create a list with several objects. The most important are:
#' - fit.y: For each MCMC sample, chi-squared statistic for total number of detections
#' for all sites in data wrt expected model expectation
#' - fit.y.rep: For each MCMC sample, chi-squared statistic for total number of detections
#' for all sites in simulated data wrt expected model expectation
#' - fit.y.group.quants: For each site, chi-squared statistic for total number of detections
#' in data wrt expected model expectation (posterior distribution quantiles)
#' - fit.y.rep.group.quants: For each site, chi-squared statistic for total number of detections
#' in simulated data wrt expected model expectation (posterior distribution quantiles)
#' - y.summ.per.site: Data frame with number of detections per site obtained from
#' data simulated from the posterior distribution. Also the observed number of detections
#' in the data.
#' @export
diagnoseGofSpOccu <- function(object, post_sims, fit_stat = "chi-squared", group = 1){

    cl <- match.call()

    # Functions
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

    out <- list()

    # Single-species models

    y <- object$y
    J <- nrow(y)

    if (is(object, 'PGOcc')) {

        ## EDITED: fitted.spPGOcc() ALREADY KILLS MY COMPUTER!!!
        # fitted.out <- simDetSpOccu(object)
        ##

    } else {
        stop("Only non-spatial models supported at the moment.")
        #fitted.out <- fitted.spPGOcc(object)
    }

    z.samples <- object$z.samples
    y.rep.samples <- post_sims$y.rep.samples
    det.prob <- post_sims$p.samples
    n.samples <- object$n.post * object$n.chains
    fit.y <- rep(NA, n.samples)
    fit.y.rep <- rep(NA, n.samples)
    e <- 0.0001

    # Do the stuff
    if (group == 1) {

        # Compute sum of detections over all visits to each site
        y.grouped <- apply(y, 1, sum, na.rm = TRUE)

        ## edited ##

        # Compute sum of detections over all visit to each site for each MCMC iteration
        # y.rep.grouped <- apply(y.rep.samples, c(1, 2), sum, na.rm = TRUE)

        y.rep.grouped <- as.data.frame(y.rep.samples) %>%
            dplyr::mutate(site_id = attr(y.rep.samples, "site_id")) %>%
            dplyr::group_by(site_id) %>%
            dplyr::summarise(dplyr::across(.cols = dplyr::starts_with("V"), .fns = ~sum(.x))) %>%
            dplyr::ungroup() %>%
            dplyr::select(-1) %>%
            as.matrix()

        # Compare number of detections observed vs simulated per pentad
        diff_dets <- sweep(y.rep.grouped, 1, y.grouped)
        mean_diff_dets <- apply(diff_dets, 1, mean)

        names(mean_diff_dets) <- names(y.grouped)



        ####

        fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
        fit.big.y <- matrix(NA, length(y.grouped), n.samples)

        ## Edited ##
        # We need to compute the latent state for each visit, rather than for each site
        z.samples.long <- t(z.samples)[attr(det.prob, "site_id"),]

        # Calculate the expected conditional detection probability per visit
        det.visit <- det.prob * z.samples.long

        # Calculate the expected conditional number of detections per site
        E.grouped <- as.data.frame(det.visit) %>%
            dplyr::mutate(site_id = attr(det.prob, "site_id")) %>%
            dplyr::group_by(site_id) %>%
            dplyr::summarise(dplyr::across(.cols = dplyr::starts_with("V"), .fns = ~sum(.x))) %>%
            dplyr::ungroup() %>%
            dplyr::select(-1) %>%
            as.matrix()


        ###

        if (fit_stat %in% c('chi-squared', 'chi-square')) {
            for (i in 1:n.samples) {

                ## Edited ##
                # For each MCMC iter simulate detection/non-detection for each visit,
                # then sum for each site
                # E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 2, sum, na.rm = TRUE)


                # For each MCMC iteration, calculate the chi-square statistic for the observed data
                fit.big.y[, i] <- (y.grouped - E.grouped[,i])^2 / (E.grouped[,i] + e)
                fit.y[i] <- sum(fit.big.y[, i])

                # For each MCMC iteration, calculate the chi-square statistic for the simulated data
                fit.big.y.rep[, i] <- (y.rep.grouped[,i] - E.grouped[,i])^2 / (E.grouped[,i] + e)
                fit.y.rep[i] <- sum(fit.big.y.rep[, i])

                bayes_p <- mean(fit.y.rep > fit.y)
                ##

            }
        } else if (fit_stat == 'freeman-tukey') {
            stop("Only chi-squared test supported at the moment")
        }
    } else if (group == 2) {
        stop("Only grouping by site supported at the moment")
    } else if (fit_stat == 'freeman-tukey') {
        stop("Only chi-squared test supported at the moment")
    }

    # Prepare output

    out$fit.y <- fit.y
    out$fit.y.rep <- fit.y.rep
    out$fit.y.group.quants <- apply(fit.big.y, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    out$fit.y.rep.group.quants <- apply(fit.big.y.rep, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))

    out$mean_diff_dets <- mean_diff_dets

    # For summaries
    out$pentad <- dimnames(object$y)[[1]]
    out$group <- group
    out$fit.stat <- fit_stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
    out$n.chains <- object$n.chains
    out$bayes_p <- round(bayes_p, 3)

    # class(out) <- 'ppcOcc'

    return(out)

}





#' Diagnose convergence for spOccupancy model
#'
#' @description This function runs basic Rhat checks for an spOccupancy model
#' @inheritParams ppl_run_pipe_abu1
#' @param fit A `spOccupancy` model fit to ABAP data.
#' @param year The year the pipeline is run for
#'
#' @return A data frame with Rhat values for the different parameters estimated
#' by the model.
#'
#' @export
#'
#' @examples
diagnoseRhatSpOccu <- function(fit, sp_code, year){

    # Extract rhat values

    # For betas
    rhats_beta <- fit$rhat$beta
    names(rhats_beta) <- paste0(colnames(fit$beta.samples), "_psi")

    # Psi random effects
    rhats_rm_psi <- fit$rhat$sigma.sq.psi
    if(!is.null(rhats_rm_psi))
        names(rhats_rm_psi) <- paste0(colnames(fit$sigma.sq.psi.samples), "_var_psi")

    # For alphas
    rhats_alpha <- fit$rhat$alpha
    names(rhats_alpha) <- paste0(colnames(fit$alpha.samples), "_p")

    # p random effects
    rhats_rm_p <- fit$rhat$sigma.sq.p
    if(!is.null(rhats_rm_p))
        names(rhats_rm_p) <- paste0(colnames(fit$sigma.sq.p.samples), "_var_p")

    rhats_vec <- c(rhats_beta, rhats_rm_psi, rhats_alpha, rhats_rm_p)
    names(rhats_vec) <- gsub("\\(|\\)", "", names(rhats_vec))

    # 0 means convergence 1 non-convergence
    rhats <- lapply(rhats_vec, function(x) sum(abs(x - 1) > 0.1))

    # Make data frame
    rhats <-  rhats %>%
        as.data.frame()

    # Number of non-convergent parameters and number of detections
    rhats$nc_pars <- sum(rhats != 0)
    rhats$ndet <- sum(apply(fit$y, 1, sum, na.rm = TRUE))

    # Total number of parameters and observations
    rhats$npars <- ncol(rhats) - 2
    rhats$nsites <- nrow(fit$X)
    rhats$nvisits <- nrow(fit$X.p)

    # General info
    rhats$sp <- sp_code
    rhats$years <- year

    return(rhats)

}




#' Fit spOccupancy model
#'
#' @param site_data_year Occupancy site data for a single year and species
#' (see \code{\link{ppl_create_site_visit}})
#' @param visit_data_year Occupancy visit data for a single year and species
#' (see \code{\link{ppl_create_site_visit}})
#' @param config A list with pipeline configuration parameters
#' (see \code{\link{configPipeline}}).
#' @param sp_code SAFRING code of the species the pipeline is running for
#' @param spatial Logical, indicating whether spatial random effects should be
#' included in the model (TRUE) or not (FALSE, default)
#' @param sp_sites Spatial object containing the pentads in `site_data_year`.
#' @param ... Other arguments that might be needed (e.g. for messages)
#'
#' @return Either a spOccupancy model fit or the integer 3, indicating that model fit
#' failed.
#' @export
#'
#' @examples
fitSpOccu <- function(site_data_year, visit_data_year, config, sp_code, spatial = FALSE, sp_sites, ...){

    # Prepare data for spOccupancy
    occu_data <- prepSpOccuData_single(site_data_year, visit_data_year, config, spatial = spatial, sp_sites)


    ## Define models

    # Priors and initial values
    # Note: we could use posterior of previous years models to define priors

    # For single season models
    # Number of samples
    n_samples <- 2e4
    batch_length <- 25
    n_batch  <- n_samples/batch_length

    if(spatial){

        # Pair-wise distances between all sites
        dist_sites <- stats::dist(occu_data$coords)/1000

        # Specify list of inits
        inits <- list(alpha = 0,
                      beta = 0,
                      z = apply(occu_data$y, 1, max, na.rm = TRUE),
                      sigma.sq = 1,
                      phi = 3 / mean(dist_sites),
                      w = rep(0, nrow(occu_data$y)))

        # Priors
        priors <- list(alpha.normal = list(mean = 0, var = 2.72),
                       beta.normal = list(mean = 0, var = 2.72),
                       sigma.sq.ig = c(2, 1),
                       phi.unif = c(3/(1000*min(dist_sites)), 3/(min(dist_sites))))

        # Run model
        fit <- spOccupancy::spPGOcc(occ.formula = reformulate(config$site_mod),
                                    det.formula = reformulate(config$visit_mod),
                                    cov.model = "exponential", NNGP = TRUE, n.neighbors = 10,
                                    data = occu_data, inits = inits, priors = priors,
                                    batch.length = batch_length, n.batch = n_batch, n.burn = 2000,
                                    accept.rate = 0.43, tuning = list(phi = 4),
                                    n.omp.threads = 3, n.thin = 20, n.chains = 3,
                                    verbose = TRUE, n.report = 200)

    } else {

        # Look for the most recent model fit within 40 years
        filenames <- setSpOutFilePath("occu_fit", config, year_sel + seq(-20, 20, 1), sp_code, ".rds")

        files_exist <- file.exists(filenames)

        if(any(files_exist)){
            filename <- filenames[max(which(files_exist))]
        } else {
            filename <- NULL
        }

        # Try default priors for now
        # filename <- NULL

        # If there is one, set priors based on those posterior estimates
        if(!is.null(filename)){

            # Load previous fit
            prev_fit <- readRDS(filename)

            # Define priors
            priors <- defineSpOccuPriors(prev_fit, config)

            # Specify list of inits
            inits <- list(alpha = priors$alpha.normal$mean,
                          beta = priors$beta.normal$mean,
                          z = apply(occu_data$y, 1, max, na.rm = TRUE))

            # Add random effects inits if necessary
            re_occ <- grepl("\\|", config$occ_mod)
            re_det <- grepl("\\|", config$det_mod)

            if(any(re_occ)){
                inits <- c(inits,
                           list(sigma.sq.psi.ig = list(shape = rep(0.1, sum(re_occ)),
                                                       scale = rep(0.1, sum(re_occ)))))
            }

            if(any(re_det)){
                inits <- c(inits,
                           list(sigma.sq.p.ig = list(shape = rep(0.1, sum(re_det)),
                                                     scale = rep(0.1, sum(re_det)))))
            }

            message(paste("Setting priors using", filename))

        } else {

            # Generic priors
            priors <- list(alpha.normal = list(mean = 0, var = 2.5),
                           beta.normal = list(mean = 0, var = 2.5))

            # Specify list of inits
            inits <- list(alpha = 0,
                          beta = 0,
                          z = apply(occu_data$y, 1, max, na.rm = TRUE))

            # Add random effects if necessary
            re_occ <- grepl("\\|", config$occ_mod)
            re_det <- grepl("\\|", config$det_mod)

            if(any(re_occ)){
                priors <- c(priors,
                            list(sigma.sq.psi.ig = list(shape = rep(1, sum(re_occ)),
                                                        scale = rep(1.5, sum(re_occ)))))
                inits <- c(inits,
                           list(sigma.sq.psi.ig = list(shape = rep(0.1, sum(re_occ)),
                                                       scale = rep(0.1, sum(re_occ)))))
            }

            if(any(re_det)){
                priors <- c(priors,
                            list(sigma.sq.p.ig = list(shape = rep(1, sum(re_det)),
                                                      scale = rep(1.5, sum(re_det)))))
                inits <- c(inits,
                           list(sigma.sq.p.ig = list(shape = rep(0.1, sum(re_det)),
                                                     scale = rep(0.1, sum(re_det)))))
            }

            message("Using generic priors")

        }


        # Run model

        fit <- tryCatch({
            out <- spOccupancy::PGOcc(occ.formula = reformulate(config$occ_mod),
                                      det.formula = reformulate(config$det_mod),
                                      data = occu_data, inits = inits, priors = priors,
                                      n.samples = n_samples, n.omp.threads = 1,
                                      n.thin = 20, n.chains = 3,
                                      verbose = TRUE, n.report = n_samples)

            out

        },
        error = function(cond) {
            sink(file.path(config$out_dir, paste0("reports/error_occu_fit_", year_sel, "_", sp_code, ".txt")))
            print(cond)
            sink()
            message(cond)
            return(NULL)
        },
        warning = function(cond) {
            sink(file.path(config$out_dir, paste0("reports/warning_occu_fit_", year_sel, "_", sp_code, ".txt")))
            print(cond)
            sink()
            message(cond)
            return(out)
        })

    }

    # Save fit and return 0 if success
    if(!is.null(fit)){

        # Save covariate scaling factors
        fit$det.scale <- list(scale = unlist(lapply(occu_data$det.covs, attr, "scaled:scale")),
                              center = unlist(lapply(occu_data$det.covs, attr, "scaled:center")))

        fit$occ.scale <- list(scale = attr(occu_data$occ.covs, "scaled:scale"),
                              center = attr(occu_data$occ.covs, "scaled:center"))

        fit$priors <- priors


        return(fit)

    } else {

        return(3)

    }

}





#' Predict from spOccupancy model fit
#'
#' @inheritParams ppl_summarise_occu
#'
#' @return A list with two elements: 1) posterior occupancy probability samples for
#' the South African pentads, and 2) posterior detection probability samples for each
#' visit in the ABAP data for the year `year_sel`.
#' @export
#'
#' @examples
predictSpOccu <- function(fit, sp_code, year_sel, config, ...){

    varargs <- list(...)

    # Prepare prediction data

    # Load site data for prediction
    sitefile <- file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv"))
    visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

    # Load data
    sitedata <- utils::read.csv(sitefile, check.names = FALSE)
    visitdata <- utils::read.csv(visitfile, check.names = FALSE)
    detdata <- utils::read.csv(detfile, check.names = FALSE)


    ## Format for occupancy modelling

    occudata <- BIRDIE::createOccuData(config = config,
                                       sp_code = sp_code,
                                       years = year_sel,
                                       site_data = sitedata,
                                       visit_data = NULL)

    # Add detection info to visit data and subset years
    occudata$visit <- visitdata %>%
        dplyr::left_join(detdata,
                         by = c("CardNo", "StartDate", "year", "Pentad")) %>%
        dplyr::filter(year == year_sel)

    # Save codes of pentads we are predicting occupancy and detection for
    occ_pentads <- occudata$site$Pentad
    det_pentads <- occudata$visit$Pentad

    # Save also detections in visits
    obs_visit <- occudata$visit$obs

    # Select occupancy variables to be included in the model
    occ_vars <- dimnames(fit$X)[[2]]
    occ_vars <- occ_vars[occ_vars != "(Intercept)"]          # Remove intercept
    occ_vars <- c(occ_vars, dimnames(fit$X.psi.re)[[2]])       # Add random effects

    occudata$site <- occudata$site %>%
        dplyr::select(dplyr::all_of(occ_vars))

    # Select detection variables to be included in the model
    det_vars <- dimnames(fit$X.p)[[2]]
    det_vars <- det_vars[det_vars != "(Intercept)"]          # Remove intercept
    det_vars <- c(det_vars, dimnames(fit$X.p.re)[[2]])       # Add random effects

    occudata$visit <- occudata$visit %>%
        dplyr::select(dplyr::all_of(det_vars))

    # Scale and center as for model
    scale_fct_occ <- fit$occ.scale
    occudata$site[,names(scale_fct_occ$center)] <- scale(occudata$site[,names(scale_fct_occ$center)],
                                                         center = scale_fct_occ$center,
                                                         scale =  scale_fct_occ$scale)

    scale_fct_det <- fit$det.scale
    occudata$visit[,names(scale_fct_det$center)] <- scale(occudata$visit[,names(scale_fct_det$center)],
                                                          center = scale_fct_det$center,
                                                          scale =  scale_fct_det$scale)

    # Add intercept
    occudata <- purrr::map(occudata, ~ .x %>%
                               dplyr::mutate(intcp = 1) %>%
                               dplyr::select(intcp, dplyr::everything()))


    pred_data <- list(psi = NA, p = NA)
    pred_data$psi <- spOccupancy:::predict.PGOcc(fit, as.matrix(occudata$site), ignore.RE = FALSE, type = "occupancy")$psi.0.samples
    pred_data$p <- spOccupancy:::predict.PGOcc(fit, as.matrix(occudata$visit), ignore.RE = FALSE, type = "detection")$p.0.samples

    # Add pentad information
    attr(pred_data$psi, "pentads") <- occ_pentads
    attr(pred_data$p, "pentads") <- det_pentads

    # And detections information
    attr(pred_data$p, "obs") <- obs_visit

    # Add year information
    attr(pred_data$psi, "year") <- year_sel
    attr(pred_data$p, "year") <- year_sel

    return(pred_data)

}



#' Prepare spOccupancy data for single-season model fitting
#'
#' @inheritParams ppl_run_pipe_dst1
#' @param site_data A data frame containing information about covariates associated
#' with ABAP pentads for a given year.
#' @param visit_data A data frame with information associated with sampling visits
#' to ABAP pentads in a given year. Detection/non-detection data for the species
#' of interest must also be included in this data frame.
#' @param spatial Logical, indicating whether spatial random effects should be
#' included in the model. Defaults to FALSE.
#' @param sp_sites Spatial object with the sites to be used for fitting spatial models
#'
#' @return
#' @export
#'
#' @examples
prepSpOccuData_single <- function(site_data, visit_data, config, spatial = FALSE, sp_sites = NULL){


    ## Prepare spOccupancy data list

    # Return list for spatial occupancy model
    if(spatial){
        # Add coordinates to the data
        spocc_data <- ABAP::abapToSpOcc_single(visit_data,
                                               pentads = sp_sites %>%
                                                   dplyr::filter(Name %in% unique(site_data$Name)))
    } else {
        spocc_data <- ABAP::abapToSpOcc_single(visit_data)
    }


    ## Add covariates to spOccupancy object

    # Select occupancy covariates and add to data list
    tt_occ <- stats::terms(stats::reformulate(config$occ_mod))

    occ_vars <- attr(tt_occ, "term.labels")
    occ_vars <- gsub(".* \\| ", "", occ_vars)

    occ_cov_sel <- site_data %>%
        dplyr::select(pentad = Pentad, dplyr::all_of(occ_vars))

    spocc_data <- spocc_data %>%
        ABAP::addEEtoSpOcc_single(ee_data = occ_cov_sel)

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "occ", scale_vars = names(occ_cov_sel)[-1])

    # Add detection covariates
    tt_det <- stats::terms(stats::reformulate(config$det_mod))

    det_vars <- attr(tt_det, "term.labels")
    det_vars <- gsub(".* \\| ", "", det_vars)

    det_cov_sel <- visit_data %>%
        dplyr::select(Pentad, StartDate, dplyr::all_of(det_vars))

    spocc_data <- addSpOccDetCovt(spocc_data, det_cov_sel)

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "det", scale_vars = c("log_hours", "prcp", "tdiff"))

    # # Create a cyclic month variable
    # spocc_data$det.covs$month_sin <- sin(2*pi*spocc_data$det.covs$month/12)
    # spocc_data$det.covs$month_cos <- cos(2*pi*spocc_data$det.covs$month/12)


    return(spocc_data)

}


#' Prepare spOccupancy data for multi-season model fitting
#'
#' @inheritParams ppl_run_pipe_dst1
#' @param spatial Logical, indicating whether spatial random effects should be
#' included in the model. Defaults to FALSE.
#'
#' @return
#' @export
#'
#' @examples
prepSpOccuData_multi <- function(sp_code, year, config, spatial = FALSE, ...){

    varargs <- list(...)


    ## Load data

    # File names
    visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
    sitefile <- file.path(config$out_dir, paste0("occu_site_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

    # Read in site and visit data
    site_data <- utils::read.csv(sitefile)
    visit_data <- utils::read.csv(visitfile)
    det_data <- utils::read.csv(detfile)

    # Stop if there are no detections
    if(!1 %in% unique(det_data$obs)){
        warning(paste("No detection of species", sp_code))
        return(1)
    }

    # Or species detected in too few Pentads
    n_pentads <- det_data %>%
        dplyr::count(Pentad, obs) %>%
        dplyr::filter(obs == 1) %>%
        nrow()

    if(n_pentads < 5){
        warning(paste("Species", sp_code, "detected in less than 5 pentads"))
        return(2)
    }


    ## Prepare spOccupancy data list

    # Add detection info to visit data
    visit_data <- visit_data %>%
        dplyr::left_join(det_data,
                         by = c("CardNo", "StartDate", "year", "Pentad")) %>%
        dplyr::mutate(Spp = ifelse(obs == 0, "-", 1))

    # Return list for spatial occupancy model
    if(spatial){
        # Add coordinates to the data
        spocc_data <- ABAP::abapToSpOcc_multi(visit_data,
                                              pentads = sa_pentads %>%
                                                  dplyr::filter(Name %in% unique(site_data$Name)),
                                              seasons = "year")
    } else {
        spocc_data <- ABAP::abapToSpOcc_multi(visit_data, seasons = "year")
    }


    ## Add covariates to spOccupancy object

    # Keep only those sites that appear in visits
    site_data <- site_data %>%
        dplyr::filter(Pentad %in% unique(visit_data$Pentad))

    # Select covariates and add to data list
    spocc_data <- spocc_data %>%
        ABAP::addEEtoSpOcc_multi(
            ee_data = site_data %>%
                dplyr::select(pentad = Pentad, year, dist_coast, prcp, tdiff, ndvi,
                              watext, watrec, elev, log_dist_coast, log_watext),
            type = "occ", seasons = "year")

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "occ", scale_vars = names(spocc_data$occ.covs)[-1])

    # Add detection covariates
    spocc_data <- spocc_data %>%
        addSpOccDetCovt(visit_data %>%
                            dplyr::mutate(StartMonth = lubridate::month(StartDate)) %>%
                            dplyr::select(Pentad, StartDate, StartMonth, year),
                        seasons = "year")


    # Add observer ID as covariate
    spocc_data <- spocc_data %>%
        addSpOccDetCovt(visit_data %>%
                            dplyr::mutate(obs_id = as.numeric(ObserverNo)) %>%
                            dplyr::select(Pentad, StartDate, obs_id, year),
                        seasons = "year")

    # Add temperature, precipitation, CWAC site presence and pentad as detection covariates
    spocc_data <- spocc_data %>%
        ABAP::addEEtoSpOcc_multi(
            ee_data = site_data %>%
                dplyr::group_by(Pentad) %>%
                dplyr::mutate(site_id = dplyr::cur_group_id()) %>%
                dplyr::ungroup() %>%
                dplyr::select(pentad = Pentad, year, prcp, tdiff, cwac, site_id),
            type = "det", seasons = "year")

    # Scale covariates
    spocc_data <- scaleSpOccVars(spocc_data, "det", scale_vars = c("prcp", "tdiff"))

    # # Create a cyclic month variable
    # spocc_data$det.covs$month_sin <- sin(2*pi*spocc_data$det.covs$month/12)
    # spocc_data$det.covs$month_cos <- cos(2*pi*spocc_data$det.covs$month/12)


    return(spocc_data)

}



#' Scale covariates in spOccupancy-type data
#'
#' @param spOcc_data an spOccupancy data list.
#' @param var_type Type of variables we want to scale. Currently, one of "occ"
#' occupancy covariates, "det" detection covariates.
#' @param scale_vars A vector with the names of the covariates that we want to
#' scale.
#'
#' @return An spOccupancy data list with the scaled covariates, substituting the
#' original, unscaled covariates. A single factor is used to center and scale all
#' data (across all dimensions) of the covariate. That means, for example, that
#' all seasons will be scaled by the same amount. The factors used
#' for centering and scaling are stored as attributes of each covariate.
#' @export
#'
#' @examples
scaleSpOccVars <- function(spOcc_data, var_type, scale_vars){

    # if(var_type == "det"){
    #     stop("Only var_type = 'occ' is supported at the moment")
    # }

    if(var_type == "occ"){

        if(is.list(spOcc_data$occ.covs)){

            f <- function(x){
                scl <- stats::sd(c(x), na.rm = TRUE)
                cnt <- mean(c(x), na.rm = TRUE)

                out_scl <- apply(x, 2, FUN = scale, center = cnt, scale = scl)
                dimnames(out_scl) <- dimnames(x)
                attr(out_scl, "scaled:scale") <- scl
                attr(out_scl, "scaled:center") <- cnt
                out_scl
            }

            covts <- spOcc_data$occ.covs[scale_vars]

            covts <- lapply(covts, f)

            spOcc_data$occ.covs[scale_vars] <- covts

        } else if(is.matrix(spOcc_data$occ.covs)){

            covt_sel <- spOcc_data$occ.covs[,scale_vars]

            scl <- apply(covt_sel, 2, sd, na.rm = TRUE)
            cnt <- apply(covt_sel, 2, mean, na.rm = TRUE)

            # make temporary object to presenve attributes
            covt_sel <- scale(covt_sel, center = cnt, scale = scl)

            spOcc_data$occ.covs[,scale_vars] <- covt_sel
            attr(spOcc_data$occ.covs, 'scaled:center') <- attr(covt_sel, 'scaled:center')
            attr(spOcc_data$occ.covs, 'scaled:scale') <- attr(covt_sel, 'scaled:scale')

        }

    } else if(var_type == "det"){

        if(is.list(spOcc_data$det.covs)){

            f <- function(x){
                scl <- stats::sd(c(x), na.rm = TRUE)
                cnt <- mean(c(x), na.rm = TRUE)

                out_scl <- apply(x, 2, FUN = scale, center = cnt, scale = scl)
                dimnames(out_scl) <- dimnames(x)
                attr(out_scl, "scaled:scale") <- scl
                attr(out_scl, "scaled:center") <- cnt
                out_scl
            }

            covts <- spOcc_data$det.covs[scale_vars]

            covts <- lapply(covts, f)

            spOcc_data$det.covs[scale_vars] <- covts

        } else if(is.matrix(spOcc_data$det.covs)){

            covt_sel <- spOcc_data$det.covs[,scale_vars]

            scl <- apply(covt_sel, 2, sd, na.rm = TRUE)
            cnt <- apply(covt_sel, 2, mean, na.rm = TRUE)

            # make temporary object to presenve attributes
            covt_sel <- scale(covt_sel, center = cnt, scale = scl)

            spOcc_data$det.covs[,scale_vars] <- covt_sel
            attr(spOcc_data$det.covs, 'scaled:center') <- attr(covt_sel, 'scaled:center')
            attr(spOcc_data$det.covs, 'scaled:scale') <- attr(covt_sel, 'scaled:scale')

        }

    }

    return(spOcc_data)

}






#' Simulate detections from spOccupancy fit
#' @description This function comes largely from the fitted.PGOcc.R from the
#' \href{https://github.com/doserjef/spOccupancy}{spOccupancy} package
#' @param object an spOccupancy fit
#'
#' @return A list with posterior detection probabilities samples and posterior detection
#' predictions samples. The results are given in a long format and as an attribute the
#' indices of the sites the samples correspond to
#' @export
simDetSpOccu <- function(object){

    # Functions
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}


    # Object
    n.post <- object$n.post * object$n.chains
    X.p <- object$X.p
    y <- object$y
    n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
    K.max <- max(n.rep)
    J <- nrow(y)
    z.long.indx <- rep(1:J, K.max)
    z.long.indx <- z.long.indx[!is.na(c(y))]
    if (nrow(X.p) == nrow(y)) {
        X.p <- do.call(rbind, replicate(ncol(y), X.p, simplify = FALSE))
        X.p <- X.p[!is.na(c(y)), , drop = FALSE]
    }
    y <- c(y)
    y <- y[!is.na(y)]
    z.samples <- object$z.samples
    alpha.samples <- object$alpha.samples

    # Create detection probability matrix (potentiall with random effects)
    tmp.samples <- matrix(0, n.post, length(y))
    if (object$pRE) {
        # Add 1 to get it to R indexing.
        X.p.re <- object$X.p.re + 1
        for (i in 1:ncol(X.p.re)) {
            tmp.samples <- tmp.samples + object$alpha.star.samples[, X.p.re[, i]]
        }
        det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) + t(tmp.samples)))
    } else {
        det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
    }

    # Simulate detections from latent states (z) and detection probabilities
    y.rep.samples <- t(apply(det.prob.samples * z.samples[, z.long.indx],
                             2, function(a) rbinom(n.post, 1, a)))

    ## edited ##

    # tmp <- array(NA, dim = c(J * K.max, n.post))
    # names.long <- which(!is.na(c(object$y)))
    # tmp[names.long, ] <- y.rep.samples
    # y.rep.samples <- array(tmp, dim = c(J, K.max, n.post))
    # y.rep.samples <- aperm(y.rep.samples, c(3, 1, 2))
    # tmp <- array(NA, dim = c(J * K.max, n.post))
    # tmp[names.long, ] <- t(det.prob.samples)
    # det.prob.samples <- array(tmp, dim = c(J, K.max, n.post))
    # det.prob.samples <- aperm(det.prob.samples, c(3, 1, 2))
    det.prob.samples <- t(det.prob.samples)

    # add site index information
    attr(y.rep.samples, "site_id") <- z.long.indx
    attr(det.prob.samples, "site_id") <- z.long.indx

    ####

    out <- list()
    out$y.rep.samples <- y.rep.samples
    out$p.samples <- det.prob.samples
    return(out)
}





#' Summarise spOccupancy occupancy estimates
#'
#' @description spOccupancy produces estimates of detection probability for each
#' visit and occupancy probabilities for each site. This functions takes these
#' predictions and creates detection probabilities, occupancy probabilities and
#' realized occupancy probabilities (occupancy probabilities conditional on
#' observed data) for each site.
#' @param pred_psi Occupancy probabilities estimated from an a spOccupancy model.
#' It must be a matrix (or mcmc object) with each row corresponding to an MCMC sample and each column
#' corresponding to a pentad. This object must contain a "pentad" attribute indicating
#' which pentad each column correspond to, and an attribute "year" indicating what
#' year the probabilities correspond to. Outputs from \code{\link{predictSpOccu}},
#' should be readily appropriate.
#' @param pred_p Detection probabilities estimated from a spOccupancy model.
#' It must be a matrix (or mcmc object) with each row corresponding to an MCMC sample and each column
#' corresponding to a visit. This object must contain a "pentad" attribute indicating
#' which pentad each column correspond to, an attribute "year" indicating what
#' year the probabilities correspond to, and an attribute "obs" indicating whether
#' the species was detected in the visit or not. Outputs from \code{\link{predictSpOccu}},
#' should be readily appropriate.
#' @param quants Quantiles to summarise predictions distribution passed as c("lower", "med", "upper").
#'
#' @return A tibble with estimates and/or quantiles for each pentad in site_data:
#' - psi: occupancy probability,
#' - p: detection probability,
#' - occu: realized occupancy (occupancy conditional on data).
#' @export
#'
#' @examples
summariseSpOccu <- function(pred_p, pred_psi, quants){

    ## Estimate realized occupancy

    # Calculate probability of non-detections for each pentad visited
    p_nondet <- data.frame(pentad = attr(pred_p, "pentads"),
                           obs = attr(pred_p, "obs")) %>%
        dplyr::mutate(ub = apply(pred_p, 2, quantile, quants[3]),
                      lb = apply(pred_p, 2, quantile, quants[1]),
                      med = apply(pred_p, 2, quantile, quants[2]),
                      est = apply(pred_p, 2, mean)) %>%
        tidyr::pivot_longer(cols = c("ub", "lb", "med", "est"),
                            names_to = "lim", values_to = "p") %>%
        dplyr::group_by(pentad, lim) %>%
        dplyr::summarise(q = prod(1-p),
                         obs = max(obs)) %>%
        dplyr::ungroup()

    # From probability of non-detection calculate the conditional occupancy probs
    # and plot
    pred_occu <- data.frame(pentad = attr(pred_psi, "pentads"),
                            year = attr(pred_psi, "year")) %>%
        dplyr::mutate(ub = apply(pred_psi, 2, quantile, quants[3]),
                      lb = apply(pred_psi, 2, quantile, quants[1]),
                      med = apply(pred_psi, 2, quantile, quants[2]),
                      est = apply(pred_psi, 2, mean)) %>%
        tidyr::pivot_longer(cols = c("ub", "lb", "med", "est"),
                            names_to = "lim", values_to = "psi") %>%
        dplyr::left_join(p_nondet %>%
                             dplyr::select(pentad, lim, q, obs), by = c("pentad", "lim")) %>%
        dplyr::mutate(real_occu = dplyr::case_when(obs == 1 ~ 1,
                                                   is.na(obs) ~ psi,
                                                   obs == 0 ~ psi*q / (1 - psi + psi*q)),
                      p = 1 - q) %>%
        dplyr::select(pentad, year, psi, p, real_occu, lim)

    return(pred_occu)

}
