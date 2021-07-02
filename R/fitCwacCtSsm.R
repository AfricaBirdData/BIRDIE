#' Fit a continuous-time state-space model to CWAC counts
#' (EXPERIMENTAL)
#'
#' @param counts A data frame with at least two columns: "counts" - an integer
#' column corresponding to the counts of a year and season and "season_id" - an
#' integer column that identifies the season (1 for summer and 2 for winter).
#' @param mod_file A character string corresponding to the directory where the
#' JAGS model lives at. If not provided, writeJagsModelFile is called to create
#' a "default" model.
#' @param param A vector with names of parameters to monitor.
#' @param jags_control A list specifying the JAGS MCMC settings. Possible
#' options are:
#' inits - initial parameter values, ni - total number of iterations,
#' nb - number of iterations to burn, nt - chain thinning, nc - number of chains,
#' na - number of adapting iterations, ncores - number of cores to use.
#'
#' @return A JAGS fit object
#' @export
#'
#' @examples
fitCwacCtSsm <- function(counts, mod_file = NULL, param, jags_control = NULL){

    if(is.null(mod_file)){
        modpath <- BIRDIE::writeJagsModelFile()
        warning("mod_file not provided, running writeJagsModelFile, see ?writeJagsModelFile")
    } else {
        modpath <- mod_file
    }

    # Prepare data
    data <- list(obs = log(counts$count + 0.1),
                 summer = counts$summer,
                 winter = counts$winter,
                 dsummer = counts$dseason[counts$summer == 1],
                 dwinter = counts$dseason[counts$winter == 1],
                 dt = counts$dt,
                 N = nrow(counts),
                 Ns = sum(counts$summer),
                 Nw = sum(counts$winter))

    # MCMC settings
    inits <- if(!is.null(jags_control$init)) jags_control$init # Inits function
    ni <- if(!is.null(jags_control$ni)) jags_control$ni else 10000 # number of iterations
    nb <- if(!is.null(jags_control$nb)) jags_control$nb else 5000 # burning iterations
    nt <- if(!is.null(jags_control$nt)) jags_control$nt else 1 # chain thinning
    nc <- if(!is.null(jags_control$nc)) jags_control$nc else 3 # number of chains
    na <- if(!is.null(jags_control$na)) jags_control$na else NULL # (default) adapting iterations
    ncores <- if(!is.null(jags_control$ncores)) jags_control$ncores else 1 # number of cores
    prll <- if(!is.null(jags_control$ncores)) TRUE else FALSE

    # Start Gibbs sampling
    fit <- jagsUI::jags(data = data, inits = inits,
                    parameters.to.save = param, model.file = modpath,
                    n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt,
                    modules = c('glm','lecuyer', 'dic'), factories = NULL, parallel = prll, n.cores = ncores,
                    DIC = TRUE, verbose = TRUE)

    # Remove temporary model file
    if(is.null(mod_file)){
        unlink(modpath)
    }

    return(fit)

}
