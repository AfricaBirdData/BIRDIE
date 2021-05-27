#' Fit a 2-season state-space model to CWAC counts
#'
#' @param counts A data frame with at least two columns: "counts" - an integer column corresponding to the counts of a year and season
#'     and "season_id" - an integer column that identifies the season (1 for summer and 2 for winter).
#' @param mod_file A character string corresponding to the directory where the JAGS model lives at.
#' @param param A vector with names of parameters to monitor.
#' @param jags_control A list specifying the JAGS MCMC settings. Possible options are: inits - initial parameter values,
#'     ni - total number of iterations, nb - number of iterations to burn, nt - chain thinning, nc - number of chains,
#'     na - number of adapting iterations, ncores - number of cores to use.
#'
#' @return A JAGS fit object
#' @export
#'
#' @examples
#' counts <- getCwacSiteCounts(26352535)
#' counts <- prepSsmData(counts)
#' fitCwacSsm2ss(counts, mod_file = "analysis/models/cwac_ssm_2ss_fxd.jags",
#' param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))
fitCwacSsm2ss <- function(counts, mod_file, param, jags_control = NULL){

    # Prepare data
    data <- list(winter = log(counts[counts$season_id == 2, "count", drop = TRUE]),
                     summer = log(counts[counts$season_id == 1, "count", drop = TRUE]),
                     N = nrow(counts)/2)

    # MCMC settings
    inits <- if(!is.null(jags_control$init)) jags_control$init # Inits function
    ni <- if(!is.null(jags_control$ni)) jags_control$ni else 10000 # number of iterations
    nb <- if(!is.null(jags_control$nb)) jags_control$nb else 5000 # burning iterations
    nt <- if(!is.null(jags_control$nt)) jags_control$nt else 1 # chain thinning
    nc <- if(!is.null(jags_control$nc)) jags_control$nb else 3 # number of chains
    na <- if(!is.null(jags_control$na)) jags_control$nb else NULL # (default) adapting iterations
    ncores <- if(!is.null(jags_control$ncores)) jags_control$ncores else 1 # number of cores
    prll <- if(!is.null(jags_control$ncores)) TRUE else FALSE

    # Start Gibbs sampling
    fit <- jags(data = data, inits = inits,
                    parameters.to.save = param, model.file = mod_file,
                    n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt,
                    modules = c('glm','lecuyer', 'dic'), factories = NULL, parallel = prll, n.cores = ncores,
                    DIC = TRUE, verbose = TRUE)

    return(fit)

}
