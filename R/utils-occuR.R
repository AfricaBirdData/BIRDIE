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
