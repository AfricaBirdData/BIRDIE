fitCwacSsm2ss <- function(counts, mod_file, jags_control = NULL){

    # Prepare data
    mod.data <- list(winter = log(counts[counts$season_id == 2, "count", drop = TRUE]),
                     summer = log(counts[counts$season_id == 1, "count", drop = TRUE]),
                     N = nrow(counts)/2)

    # MCMC settings
    mod.inits <- if(!is.null(jags_control$init)) jags_control$init # Inits function
    ni <- if(!is.null(jags_control$ni)) jags_control$ni else 10000 # number of iterations
    nb <- if(!is.null(jags_control$nb)) jags_control$nb else 5000 # burning iterations
    nt <- if(!is.null(jags_control$nt)) jags_control$nt else 1 # chain thinning
    nc <- if(!is.null(jags_control$nc)) jags_control$nb else 3 # number of chains
    na <- if(!is.null(jags_control$na)) jags_control$nb else NULL # (default) adapting iterations
    ncores <- if(!is.null(jags_control$ncores)) jags_control$ncores else 1 # number of cores
    prll <- if(!is.null(jags_control$ncores)) TRUE else FALSE

    # Parameters to estimate
    mod.param <- c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e")
    mod.param<- "mu_t"

    # Start Gibbs sampling
    mod.fit <- jags(data = mod.data, inits = mod.inits,
                    parameters.to.save = mod.param, model.file = "analysis/models/cwac_ssm_2season.jags",
                    n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt,
                    modules = c('glm','lecuyer', 'dic'), factories = NULL, parallel = prll, n.cores = ncores,
                    DIC = TRUE, verbose = TRUE)

    return(mod.fit)

}
