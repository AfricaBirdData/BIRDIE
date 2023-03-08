#' Simulate detections from spOccupancy fit
#' @description This function comes largely from the fitted.PGOcc.R from the
#' \href{https://github.com/doserjef/spOccupancy}{spOccupancy} package
#' @param object an spOccupancy fit
#'
#' @return A list with posterior detection probabilities samples and posterior detection
#' predictions samples. The results are given in a long format and as an attribute the
#' indices of the sites the samples correspond to
#' @export
#' @keywords internal
#' @noRd
simDetSpOccu <- function(object){

    # Functions -------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}


    # Object ----------------------------
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
#' data. One statistic for is computed for each MCMC iteration
#' @export
#' @keywords internal
#' @noRd
gofSpOccupancy <- function(object, post_sims, fit_stat = "chi-squared", group = 1){

    cl <- match.call()

    # Functions -------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

    out <- list()

    # Single-species models -------------------------------------------------

    y <- object$y
    J <- nrow(y)
    if (is(object, 'PGOcc')) {

        ## EDITED: fitted.spPGOcc() ALREADY KILLS MY COMPUTER!!!
        # fitted.out <- simDetSpOccu(object)
        #####

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
out$n.burn <- object$n.burn
out$n.thin <- object$n.thin
out$n.post <- object$n.post
out$n.chains <- object$n.chains

# class(out) <- 'ppcOcc'

return(out)

}


defineSpOccupancyPriors <- function(prev_fit){

    prior_mean_beta <- apply(prev_fit$beta.samples, 2, mean)
    prior_sd_beta <- apply(prev_fit$beta.samples, 2, sd) + 1
    names(prior_mean_beta) <- names(prior_sd_beta) <- dimnames(prev_fit$beta.samples)[[2]]

    prior_mean_alpha <- apply(prev_fit$alpha.samples, 2, mean)
    prior_sd_alpha <- apply(prev_fit$alpha.samples, 2, sd) + 1
    names(prior_mean_alpha) <- names(prior_sd_alpha) <- dimnames(prev_fit$alpha.samples)[[2]]

    if(!is.null(prev_fit$sigma.sq.p.samples)){
        prior_mean_sigma.sq.p <- apply(prev_fit$sigma.sq.p.samples, 2, mean)
        prior_sd_sigma.sq.p <- apply(prev_fit$sigma.sq.p.samples, 2, sd) + 1
        names(prior_mean_sigma.sq.p) <- names(prior_sd_sigma.sq.p) <- dimnames(prev_fit$sigma.sq.p.samples)[[2]]

        scale.p.ig <- prior_sd_sigma.sq.p^2 / prior_mean_sigma.sq.p
        shape.p.ig <- prior_mean_sigma.sq.p / scale.p.ig

    } else {
        scale.p.ig <- NULL
        shape.p.ig <- NULL
    }

    if(!is.null(prev_fit$sigma.sq.psi.samples)){
        prior_mean_sigma.sq.psi <- apply(prev_fit$sigma.sq.psi.samples, 2, mean)
        prior_sd_sigma.sq.psi <- apply(prev_fit$sigma.sq.psi.samples, 2, sd) + 1
        names(prior_mean_sigma.sq.psi) <- names(prior_sd_sigma.sq.psi) <- dimnames(prev_fit$sigma.sq.psi.samples)[[2]]

        scale.psi.ig <- prior_sd_sigma.sq.psi^2 / prior_mean_sigma.sq.psi
        shape.psi.ig <- prior_mean_sigma.sq.psi / scale.psi.ig

    } else {
        scale.psi.ig <- NULL
        shape.psi.ig <- NULL
    }

    priors <- list(alpha.normal = list(mean = prior_mean_alpha,
                                       var = prior_sd_alpha^2),
                   beta.normal = list(mean = prior_mean_beta,
                                      var = prior_sd_beta^2))

    # Add random effects
    if(is.null(shape.psi.ig)){
        priors <- c(priors,
                    list(sigma.sq.psi.ig = list(shape = shape.psi.ig,
                                                scale = scale.psi.ig)))
    }

    if(is.null(shape.p.ig)){
        priors <- c(priors,
                    list(sigma.sq.p.ig = list(shape = shape.p.ig,
                                              scale = scale.p.ig)))
    }

    return(priors)

}

