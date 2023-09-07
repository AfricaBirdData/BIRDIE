#' Simulate detections from occuR fit
#' @description This function predicts probabilities of occupancy and detection
#' from an occuR model. From those it simulates detection/non-detection data
#' @param obj an occuR fit
#'
#' @return A list with posterior detection probabilities samples and posterior detection
#' predictions samples. The results are given in a long format and as an attribute the
#' indices of the sites the samples correspond to
#' @export
#' @keywords internal
#' @noRd
simDetOccuR <- function(obj){

    fix <- obj$res$par.fixed
    ran <- obj$res$par.random
    if(obj$mats$U_psi_n[1] > 0) obj$mats$U_psi <- as(obj$mats$U_psi, "matrix")
    if(obj$mats$U_p_n[1] > 0) obj$mats$U_p <- as(obj$mats$U_p, "matrix")
    pred <- occuR:::get_predicted_values(fix, ran, obj$mats)

    nboot <- 1000
    Q <- obj$res$jointPrecision
    if (!is.null(Q)) {
        Q <- Q[!grepl("log_lambda_|lsig_gamma_", colnames(Q)),
               !grepl("log_lambda_|lsig_gamma_", colnames(Q)), drop = FALSE]
        V <- Matrix::solve(Q)
        # Sometimes V doesn't pass the symmetry test
        if(isSymmetric(as.matrix(V), tol = 1e-4)){
            V <- Matrix::forceSymmetric(V)
        }
    } else {
        V <- obj$res$cov.fixed
    }
    param <- c(fix, ran)
    param <- param[!grepl("log_lambda|lsig_gamma_", names(param))]
    boots <- mgcv::rmvn(nboot, param, as.matrix(V))
    colnames(boots) <- names(param)
    nfix <- length(fix[!grepl("log_lambda|lsig_gamma_", names(fix))])
    boots_psi <- matrix(0, nr = nboot, nc = length(pred$psi))
    boots_p <- matrix(0, nr = nboot, nc = length(pred$p))
    boots_z <- matrix(0, nr = nboot, nc = length(pred$psi))
    boots_det <- matrix(0, nr = nboot, nc = length(pred$p))
    for (b in 1:nboot) {
        boot_res <- occuR:::get_predicted_values(boots[b, 1:nfix], boots[b, -(1:nfix)], obj$mats)
        boots_psi[b,] <- boot_res$psi
        boots_p[b,] <- boot_res$p
        boots_z[b,] <- rbinom(ncol(boots_z), 1, boot_res$psi)
    }

    # pred$psiboot <- boots_psi
    # pred$pboot <- boots_p

    # Functions -------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}


    # # Object ----------------------------
    # n.post <- obj$n.post * obj$n.chains
    # X.p <- obj$X.p
    # y <- obj$y
    # n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
    # K.max <- max(n.rep)
    # J <- nrow(y)
    # z.long.indx <- rep(1:J, K.max)
    # z.long.indx <- z.long.indx[!is.na(c(y))]
    # if (nrow(X.p) == nrow(y)) {
    #     X.p <- do.call(rbind, replicate(ncol(y), X.p, simplify = FALSE))
    #     X.p <- X.p[!is.na(c(y)), , drop = FALSE]
    # }
    # y <- c(y)
    # y <- y[!is.na(y)]
    # z.samples <- obj$z.samples
    # alpha.samples <- obj$alpha.samples
    #
    # # Create detection probability matrix (potentiall with random effects)
    # tmp.samples <- matrix(0, n.post, length(y))
    # if (obj$pRE) {
    #     # Add 1 to get it to R indexing.
    #     X.p.re <- obj$X.p.re + 1
    #     for (i in 1:ncol(X.p.re)) {
    #         tmp.samples <- tmp.samples + obj$alpha.star.samples[, X.p.re[, i]]
    #     }
    #     det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) + t(tmp.samples)))
    # } else {
    #     det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
    # }

    # Simulate detections from latent states (z) and detection probabilities
    boots_z <- boots_z[,attr(obj$mats$X_psi, "site")]
    y.rep.samples <- t(apply(boots_p * boots_z[,attr(obj$mats$X_p, "site")],
                             2, function(a) rbinom(nboot, 1, a)))

    # y.rep.samples <- t(apply(det.prob.samples * z.samples[, z.long.indx],
    #                          2, function(a) rbinom(n.post, 1, a)))

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
    boots_p <- t(boots_p)
    boots_psi <- t(boots_psi)

    # add site index information
    attr(y.rep.samples, "site_id") <- attr(obj$mats$X_p, "site")
    attr(boots_p, "site_id") <- attr(obj$mats$X_p, "site")
    attr(boots_psi, "site_id") <- attr(obj$mats$X_psi, "site")

    ####

    out <- list()
    out$y.rep.samples <- y.rep.samples
    out$p.samples <- boots_p
    out$psi.samples <- boots_psi
    out$z.samples <- t(boots_z)
    return(out)
}


#' Conduct goodness-of-fit test for occuR occupancy models
#' @description This function comes largely from the ppcOcc.R from the
#' \href{https://github.com/doserjef/spOccupancy}{spOccupancy} package
#' @param object an occuR fit
#' @param post_sims A list of posterior simulations .
#' @param fit_stat Goodness of fit statistic to compute. Currently only Chi-squared
#' is supported.
#' @param group Whether the GOF statistic should be computed for sites (1) or
#' for visits (2). Currently, only grouping by site is supported
#'
#' @return A list GOF statistics computed for observed data and for simulated
#' data. One statistic for is computed for each MCMC iteration
#' @export
#' @keywords internal
#' @noRd
gofOccuR <- function(object, post_sims, fit_stat = "chi-squared", group = 1){

    cl <- match.call()

    # Functions -------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

    out <- list()

    # Single-species models -------------------------------------------------

    y <- object$mats$y
    J <- length(y)

    z.samples <- post_sims$z.samples
    y.rep.samples <- post_sims$y.rep.samples
    det.prob <- post_sims$p.samples
    occ.prob <- post_sims$psi.samples
    n.samples <- ncol(det.prob)
    fit.y <- rep(NA, n.samples)
    fit.y.rep <- rep(NA, n.samples)
    e <- 0.0001

    # Do the stuff
    if (group == 1) {

        # Compute sum of detections over all visits to each site
        y.grouped <- as.data.frame(y) %>%
            dplyr::mutate(site_id = attr(object$mats$X_p, "site")) %>%
            dplyr::group_by(site_id) %>%
            dplyr::summarise(n = sum(y)) %>%
            dplyr::ungroup() %>%
            dplyr::select(-1) %>%
            as.matrix()

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
        ####

        fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
        fit.big.y <- matrix(NA, length(y.grouped), n.samples)

        ## Edited ##
        # We need to compute the latent state for each visit, rather than for each site
        z.samples.long <- z.samples[attr(det.prob, "site_id"),]

        # Calculate the expected conditional detection probability per visit
        # det.visit <- det.prob * z.samples.long
        det.visit <- det.prob*occ.prob[attr(det.prob, "site_id"),]

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
#
#
#                 sq1 <- (y.grouped[,1] - E.grouped[,i])^2
#                 sq2 <- (y.rep.grouped[,i] - E.grouped[,i])^2

                # For each MCMC iteration, calculate the chi-square statistic for the observed data
                fit.big.y[, i] <- (y.grouped[,1] - E.grouped[,i])^2 / (E.grouped[,i] + e)
                fit.y[i] <- sum(fit.big.y[, i])

                # For each MCMC iteration, calculate the chi-square statistic for the simulated data
                fit.big.y.rep[, i] <- (y.rep.grouped[,i] - E.grouped[,i])^2 / (E.grouped[,i] + e)
                fit.y.rep[i] <- sum(fit.big.y.rep[, i])
                #####

            }
        } else if (fit_stat == 'freeman-tukey') {
            stop("Only chi-squared test supported at the moment")
        }
    } else if (group == 2) {
        stop("Only grouping by site supported at the moment")
    } else if (fit_stat == 'freeman-tukey') {
        stop("Only chi-squared test supported at the moment")
    }

out$fit.y <- fit.y
out$fit.y.rep <- fit.y.rep
out$fit.y.group.quants <- apply(fit.big.y, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
out$fit.y.rep.group.quants <- apply(fit.big.y.rep, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))

# For summaries
out$site_id <- dimnames(object$y)[[1]]
out$group <- group
out$fit.stat <- fit_stat
out$class <- class(object)
out$call <- cl
out$n.samples <- object$n.samples

# class(out) <- 'ppcOcc'

return(out)

}


#' Run diagnostics for occuR model
#'
#' @description This function runs basic diagnostic checks for occuR models.
#' Routines are based on those in spOccupancy and perform posterior predictive
#' goodness of fit checks see \code{\link[spOccupancy]{ppcOcc}}. For occuR predictions
#' are obtained via bootstrap.
#' @inheritParams ppl_run_pipe_dst1
#' @param fit An occuR model fit see \code{\link[occuR]{fit_occu}}
#' @param year_sel The year the model was run for.
#'
#' @return The function will create a list with several objects. The most important are:
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
#'
#' @export
#'
#' @examples
diagnoseOccuR <- function(fit, sp_code, config, year_sel){

    # Check convergence
    fails <- fit$fit$convergence

    if(fails != 0){
        conv_file <- file.path(config$out_dir, paste0("reports/no_converge_occu_fit_", config$package, "_", year_sel, "_", sp_code, ".txt"))
        sink(conv_file)
        print(paste("non-convergence", fit$fit$message))
        sink()
        message(paste("non-convergence", fit$fit$message)) # to console
    }

    # Simulate detections/non-detection data from model
    post_sims <- simDetOccuR(fit)

    # Check goodness-of-fit (Bayesian p-value from posterior predictive check)
    ppc_out <- gofOccuR(fit, post_sims, fit_stat = 'chi-squared', group = 1)

    bayes_p <- mean(ppc_out$fit.y.rep > ppc_out$fit.y)

    if(bayes_p < 0.05){
        sink(file.path(config$out_dir, paste0("reports/gof_occu_fit_", config$package, "_", year_sel, "_", sp_code, ".txt")))
        print(paste("Bayesian p-value =", round(bayes_p, 2)))
        sink()
        message(paste("Bayesian p-value =", round(bayes_p, 2)))
    }

    # Posterior coverage
    y_sims_per_site <- as.data.frame(post_sims$y.rep.samples) %>%
        dplyr::mutate(site_id = attr(post_sims$y.rep.samples, "site_id")) %>%
        dplyr::group_by(site_id) %>%
        dplyr::summarise(dplyr::across(.cols = dplyr::starts_with("V"), .fns = ~sum(.x))) %>%
        dplyr::ungroup() %>%
        dplyr::select(-1) %>%
        as.matrix()

    y_sims_per_site_quants <- apply(y_sims_per_site, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))

    y_grouped <- data.frame(y = fit$mats$y) %>%
        dplyr::mutate(site_id = attr(fit$mats$X_p, "site")) %>%
        dplyr::group_by(site_id) %>%
        dplyr::summarise(n = sum(y)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-1) %>%
        as.matrix()

    pst_df <- as.data.frame(t(y_sims_per_site_quants)) %>%
        stats::setNames(paste0("q", c(0.025, 0.25, 0.5, 0.75, 0.975)*100)) %>%
        dplyr::mutate(obs = y_grouped,
                      site_id = names(y_grouped))

    # Save per site info
    ppc_out$y.summ.per.site <- pst_df

    return(ppc_out)

}

#' Summarise occuR occupancy estimates
#'
#' @description occuR produces estimates of detection probability for each
#' visit and occupancy probabilities for each site. This functions takes these
#' predictions and creates detection probabilities, occupancy probabilities and
#' realized occupancy probabilities (occupancy probabilities conditional on
#' observed data) for each site.
#' @param pred_psi Occupancy probabilities estimated from an a occuR model.
#' It must be a matrix with each row corresponding to an bootstrap sample and each column
#' corresponding to a pentad. This object must contain a "pentad" attribute indicating
#' which pentad each column correspond to, and an attribute "year" indicating what
#' year the probabilities correspond to. Outputs from \code{\link{predictOccuR}},
#' should be readily appropriate.
#' @param pred_p Detection probabilities estimated from a occuR model.
#' It must be a matrix with each row corresponding to an bootstrap sample and each column
#' corresponding to a visit. This object must contain a "pentad" attribute indicating
#' which pentad each column correspond to, an attribute "year" indicating what
#' year the probabilities correspond to, and an attribute "obs" indicating whether
#' the species was detected in the visit or not. Outputs from \code{\link{predictOccuR}},
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
summariseOccuR <- function(pred_p, pred_psi, quants){

    # Estimate realized occupancy ---------------------------------------------

    # Calculate probability of non-detections for each pentad visited
    p_nondet <- data.frame(pentad = attr(pred_p, "pentad"),
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
    pred_occu <- data.frame(pentad = attr(pred_psi, "pentad"),
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



#' Predict from occuR model fit
#'
#' @inheritParams ppl_summarise_occu
#'
#' @return A list with two elements: 1) posterior occupancy probability samples for
#' the South African pentads, and 2) posterior detection probability samples for each
#' visit in the ABAP data for the year `year_sel`.
#'
#' @return
#' @export
#'
#' @examples
predictOccuR <- function(fit, sp_code, year, config, ...){

    varargs <- list(...)

    # Prepare prediction data -------------------------------------------------

    # Load site data for prediction
    sitefile <- file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv"))
    visitfile <- file.path(config$out_dir, paste0("occu_visit_dat_sa_", config$years_ch, ".csv"))
    detfile <- file.path(config$out_dir, sp_code, paste0("occu_det_dat_sa_", config$years_ch, ".csv"))

    # Load data
    sitedata <- utils::read.csv(sitefile, check.names = FALSE)
    visitdata <- utils::read.csv(visitfile, check.names = FALSE)
    detdata <- utils::read.csv(detfile, check.names = FALSE)


    # Format for occupancy modelling ------------------------------------------

    occudata <- BIRDIE::createOccuData(sp_code = sp_code,
                                       years = year_sel,
                                       site_data = sitedata,
                                       visit_data = NULL,
                                       config = config,
                                       force_abap_dwld = FALSE)

    # Add detection info to visit data and subset years
    occudata$visit <- visitdata %>%
        dplyr::left_join(detdata,
                         by = c("CardNo", "StartDate", "year", "Pentad")) %>%
        dplyr::filter(year == year_sel) %>%
        dplyr::mutate(Spp = ifelse(obs == 0, "-", "1"))

    # # Save codes of pentads we are predicting occupancy and detection for
    # occ_pentads <- occudata$site$Pentad
    # det_pentads <- occudata$visit$Pentad
    #
    # # Save also detections in visits
    # obs_visit <- occudata$visit$obs

    # Format to occuR
    occuRdata <- prepOccuRData(occudata$site, occudata$visit, config, spatial = FALSE,
                               sp_sites = NULL, scale = FALSE, keep_sites = TRUE)

    # # Select occupancy variables to be included in the model
    # tt_occ <- stats::terms(stats::reformulate(config$occ_mod))
    # occ_vars <- attr(tt_occ, "term.labels")
    # occ_vars <- gsub(".* \\| ", "", occ_vars)
    #
    # occuRdata$site <- occuRdata$site %>%
    #     dplyr::select(dplyr::all_of(occ_vars))
    #
    # # Select detection variables to be included in the model
    # tt_det <- stats::terms(stats::reformulate(config$det_mod))
    # det_vars <- attr(tt_det, "term.labels")
    # det_vars <- gsub(".* \\| ", "", det_vars)
    #
    # occuRdata$visit <- occuRdata$visit %>%
    #     dplyr::select(dplyr::all_of(det_vars))

    # Scale and center as for model
    occuRdata <- lapply(occuRdata, as.data.frame)
    scale_fct_occ <- fit$occ.scale
    scale_fct_occ <- lapply(scale_fct_occ, function(x) x[-grep(":", names(x))])   # For occuR we need to remove interactions
    occuRdata$site[,names(scale_fct_occ$center)] <- scale(occuRdata$site[,names(scale_fct_occ$center)],
                                                          center = scale_fct_occ$center,
                                                          scale =  scale_fct_occ$scale)

    scale_fct_det <- fit$det.scale
    scale_fct_det <- lapply(scale_fct_det, function(x) x[-grep(":", names(x))])   # For occuR we need to remove interactions
    occuRdata$visit[,names(scale_fct_det$center)] <- scale(occuRdata$visit[,names(scale_fct_det$center)],
                                                           center = scale_fct_det$center,
                                                           scale =  scale_fct_det$scale)

    occuRdata <- lapply(occuRdata, data.table::as.data.table)

    occuRdata$site <- occuRdata$site %>%
        dplyr::mutate(site_id = factor(site_id))
    occuRdata$visit <- occuRdata$visit %>%
        dplyr::mutate(obs_id = factor(obs_id),
                      site_id = factor(site_id))


    tryCatch({

        # Predict
        pred_occu <- predict(fit, occuRdata$visit, occuRdata$site,
                             include_re = TRUE, new_levels = TRUE, nboot = 1000)

    }, error = function(e){
        sink(file.path(config$out_dir, sp_code, paste0("failed_pred_", config$package, "_", sp_code,"_", year, ".txt")))
        print(e)
        sink()}) # TryCatch predict

    # Add pentad information
    # attr(pred_data$psi, "pentads") <- occ_pentads
    # attr(pred_data$p, "pentads") <- det_pentads

    # And detections information
    # attr(pred_data$p, "obs") <- obs_visit

    # Add year information
    # attr(pred_data$psi, "year") <- year_sel
    # attr(pred_data$p, "year") <- year_sel

    attr(pred_occu$psi, "pentad") <- occuRdata$site$pentad
    attr(pred_occu$psi, "year") <- occuRdata$site$year
    attr(pred_occu$p, "pentad") <- occuRdata$visit$pentad
    attr(pred_occu$p, "obs") <- occuRdata$visit$obs
    attr(pred_occu$psiboot, "pentad") <- occuRdata$site$pentad
    attr(pred_occu$psiboot, "year") <- occuRdata$site$year
    attr(pred_occu$pboot, "pentad") <- occuRdata$visit$pentad
    attr(pred_occu$pboot, "obs") <- occuRdata$visit$obs

    return(pred = pred_occu)

}



#' Prepare occuR data for model fitting
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
#' @param scale If TRUE covariates will be scaled. Scale factors will be provided
#' as arguments to output. Defaults to FALSE.
#' @param keep_sites If TRUE (default), all sites in `site_data` will be retained in the output.
#' This is useful for predicting from a fitted model. If FALSE, then only those
#' sites present in `visit_data` are retained, which is useful for fitting models.
#'
#' @return
#' @export
#'
#' @examples
prepOccuRData <- function(site_data, visit_data, config, spatial = FALSE,
                          sp_sites = NULL, scale = FALSE, keep_sites = TRUE){


    # Prepare occuR data list -------------------------------------------------

    # Return list for spatial occupancy model
    if(spatial){
        # Add coordinates to the data
        occur_data <- ABAP::abapToOccuR(visit_data,
                                        occasion = "year",
                                        pentads = sp_sites %>%
                                            dplyr::filter(Name %in% unique(visit_data$Name)))
    } else {
        occur_data <- ABAP::abapToOccuR(visit_data,
                                        occasion = "year")
    }


    # Add covariates to occuR object ------------------------------------------

    # Select occupancy covariates and add to data list
    tt_occ <- stats::terms(stats::reformulate(config$occ_mod))

    occ_vars <- attr(tt_occ, "term.labels")
    occ_vars <- gsub(".* \\| ", "", occ_vars)

    # Add site_id if not in covariates
    if(!"site_id" %in% occ_vars){occ_vars <- c("site_id", occ_vars)}

    occ_cov_sel <- site_data %>%
        dplyr::select(pentad = Pentad, year, dplyr::all_of(occ_vars))

    if(keep_sites){
        occur_data$site <- occ_cov_sel %>%
            dplyr::left_join(occur_data$site %>%
                                 dplyr::select(pentad, site) %>%
                                 dplyr::distinct(), by = "pentad") %>%
            dplyr::left_join(occur_data$site %>%
                                 dplyr::select(year, occasion) %>%
                                 dplyr::distinct(), by = "year") %>%
            data.table::as.data.table()
    } else {
        occur_data$site <- occur_data$site %>%
            dplyr::left_join(occ_cov_sel, by = c("pentad", "year"))
    }


    if(scale){
        # Scale covariates
        occur_data$site <- occur_data$site %>%
            dplyr::mutate(dplyr::across(-(where(is.integer) | where(is.character)),
                                        ~ scale(.x)))
    }

    # Add detection covariates
    tt_det <- stats::terms(stats::reformulate(config$det_mod))

    det_vars <- attr(tt_det, "term.labels")
    det_vars <- gsub(".* \\| ", "", det_vars)

    # Add site_id if not in covariates
    if(!"site_id" %in% det_vars){det_vars <- c("site_id", det_vars)}

    det_cov_sel <- visit_data %>%
        dplyr::arrange(Pentad, StartDate) %>%  # This is how abapToOccuR orders rows
        dplyr::select(StartDate, dplyr::all_of(det_vars))

    if(nrow(det_cov_sel) == nrow(occur_data$visit)){
        occur_data$visit <- occur_data$visit %>%
            dplyr::select(pentad, year, site, occasion, visit, obs) %>%
            dplyr::bind_cols(det_cov_sel)
    } else {
        stop("Occupancy visit data and covariate data have different number of rows")
    }

    if(scale){
        # Scale covariates
        occur_data$visit <- occur_data$visit %>%
            dplyr::mutate(dplyr::across(-(where(is.integer) | where(is.character)),
                                        ~ scale(.x)))
    }


    # # Create a cyclic month variable
    # spocc_data$det.covs$month_sin <- sin(2*pi*spocc_data$det.covs$month/12)
    # spocc_data$det.covs$month_cos <- cos(2*pi*spocc_data$det.covs$month/12)


    return(occur_data)

}



#' Select occuR model
#'
#' @param forms A list with formulas to fit an occuR model see \link[occuR]{fit_occu}
#' @param type One of "visit" if the model selection is for the probability of
#' detection or "site" if model selection is for occupancy probability. An
#' overall selection to come...
#' @param visit_data Occupancy visit data produced by \link{prepDataOccuR} or with
#' the format prescribed by \link[occuR]{fit_occu}.
#' @param site_data Occupancy site data produced by \link{prepDataOccuR} or with
#' the format prescribed by \link[occuR]{fit_occu}.
#' @param mod_id Model identifier. It is recommended to include something that
#' allows us to identify which model and to what species.
#' species id.
#'
#' @return A dataframe with four columns: df - estimated degrees of freedom,
#' form - formula of the assessed model, AIC - AIC score of the model, mod -
#' model identifier.
#' @export
#'
#' @examples
selOccuRmod <- function(forms, type, visit_data, site_data, mod_id){

    if(type == "visit"){
        form <- forms[[1]]
    } else if(type == "site"){
        form <- forms[[2]]
    }

    fit <- fit_occu(forms = forms,
                    visit_data = visit_data,
                    site_data = site_data)

    out <- data.frame(df = dof.occuR(fit, each = FALSE),
                      form = Reduce(paste, deparse(form)),
                      AIC = AIC(fit),
                      mod = mod_id)

    return(out)

}


#' Fit occuR occupancy model
#'
#' @param site_data_year Occupancy site data for a single year and species
#' (see \code{\link{ppl_create_site_visit}})
#' @param visit_data_year Occupancy visit data for a single year and species
#' (see \code{\link{ppl_create_site_visit}})
#' @param config A list with pipeline configuration parameters
#' (see \code{\link{configPreambOccu}}).
#' @param spatial Logical, indicating whether spatial random effects should be
#' included in the model (TRUE) or not (FALSE, default)
#' @param sp_sites Spatial object containing the pentads in `site_data_year`.
#' @param verbose Logical indicating whether information should be printed
#' during model fitting.
#'
#' @return Either an occuR model fit or the integer 3, indicating that model fit
#' failed.
#' @export
#'
#' @examples
fitOccuR <- function(site_data_year, visit_data_year, config, spatial = FALSE, sp_sites, verbose){

    # Prepare data for occuR
    occu_data <- prepOccuRData(site_data_year, visit_data_year, config, spatial = spatial,
                               sp_sites, scale = TRUE, keep_sites = FALSE)


    # Define models -----------------------------------------------------------

    # forms, visit_data, site_data, start = NULL, print = TRUE
    occu_data$site <- occu_data$site %>%
        dplyr::mutate(site_id = factor(site_id))
    occu_data$visit <- occu_data$visit %>%
        dplyr::mutate(site_id = factor(site_id),
                      obs_id = factor(obs_id))

    # Run model
    tryCatch({
        fit <- occuR::fit_occu(forms = list(stats::reformulate(config$occ_mod, response = "psi"),
                                            stats::reformulate(config$det_mod, response = "p")),
                               visit_data = occu_data$visit,
                               site_data = occu_data$site,
                               print = verbose)

        if(fit$fit$convergence != 0){
            warning(paste(config$package, "model failed once to converge for species", sp_code, "and year", year_sel))
            fp <- fit$res$par.fixed
            rp <- fit$res$par.random
            start_vals <- list(beta_psi = fp[names(fp) == "beta_psi"],
                               beta_p = fp[names(fp) == "beta_p"],
                               z_psi = rp[names(rp) == "z_psi"],
                               z_p = rp[names(rp) == "z_p"],
                               log_lambda_psi = fp[names(fp) == "log_lambda_psi"],
                               log_lambda_p = fp[names(fp) == "log_lambda_p"],
                               gamma_psi = rp[names(rp) == "gamma_psi"],
                               gamma_p = rp[names(rp) == "gamma_p"],
                               lsig_gamma_psi = fp[names(fp) == "lsig_gamma_psi"],
                               lsig_gamma_p = fp[names(fp) == "lsig_gamma_p"])

            start_vals[sapply(start_vals, function(x) length(x) == 0)] <- 0

            fit <- occuR::fit_occu(forms = list(stats::reformulate(config$occ_mod, response = "psi"),
                                                stats::reformulate(config$det_mod, response = "p")),
                                   visit_data = occu_data$visit,
                                   site_data = occu_data$site,
                                   start = start_vals,
                                   print = verbose)

        }


        # if(fit$fit$convergence != 0){
        #     warning(paste(config$package, "model failed twice to converge for species", sp_code, "and year", year_sel))
        #     fit <- occuR::fit_occu(forms = list(stats::reformulate(config$occ_mod, response = "psi"),
        #                                         stats::reformulate(config$det_mod[-2], response = "p")),
        #                            visit_data = occu_data$visit,
        #                            site_data = occu_data$site,
        #                            print = verbose)
        #
        # }

    },
    error = function(cond) {
        filename <- paste0("reports/error_occu_fit_", config$package, "_", year_sel, "_", sp_code, ".txt")
        sink(file.path(config$out_dir, filename))
        print(cond)
        sink()
        message(cond)
    })

    # Save fit and return 0 if success
    if(!is.null(fit)){

        # Save covariate scaling factors
        fit$det.scale <- list(scale = unlist(lapply(occu_data$visit, attr, "scaled:scale")),
                              center = unlist(lapply(occu_data$visit, attr, "scaled:center")))

        fit$occ.scale <- list(scale = unlist(lapply(occu_data$site, attr, "scaled:scale")),
                              center = unlist(lapply(occu_data$site, attr, "scaled:center")))


        return(fit)

    } else {

        return(3)

    }

}
