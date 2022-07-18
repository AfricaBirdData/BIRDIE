# This model considers two potentially unobserved states: abundance in summer and winter
# It estimates these states for a number of sites, using counts taken at those sites
# Some of the sites may be missing counts for some of the years.
# State changes are estimated in a hierarchical way such that changes across sites
# in a given year come from a common distribution.

model {

    ### PRIORS: ###

    for(s in 1:nsites){

        # Priors for initial states
        stt_s[s, 1] ~ dnorm(0, 1/10)       # prior for initial log population size
        beta[s, 1] ~ dnorm(0, 1/10)        # prior for initial log growth rate
        lambda[s, 1] ~ dnorm(0, 1/10)      # prior for log summer to winter ratio

        tau.alpha[s] ~ dgamma(2, 2)
        tau.e[s] ~ dgamma(2, 2)

        sig.alpha[s] = sqrt(1/tau.alpha[s])   # observer error standard deviation summer
        sig.e[s] = sqrt(1/tau.e[s])           # observer error standard deviation winter

        # Prior for season covariate coefficients
        for(k in 1:K){
            B[s, k] ~ dnorm(0, 1)
        }

    }

    ### latent states ###

    # State updates
    for(s in 1:nsites){
        stt_w[s, 1] = stt_s[s, 1] + lambda[s, 1]
    }

    for (y in 1:(nyears-1)){

        # tau.zeta[y] ~ dscaled.gamma(5, 2)    # same as Student-t with 2 degrees of freedom and standard deviation 5 for sig.zeta!
        tau.zeta[y] ~ dgamma(1.8, 1)
        sig.zeta[y] = sqrt(1/tau.zeta[y])

        tau.eps[y] ~ dscaled.gamma(5, 2)     # same as Student-t with 2 degrees of freedom and standard deviation 5 for sig.eps!
        # tau.eps[y] ~ dgamma(1.8, 1)
        sig.eps[y] = sqrt(1/tau.eps[y])

        for(s in 1:nsites){

            # Sample zeta and eps
            zeta[s, y] ~ dnorm(0, tau.zeta[y])
            eps[s, y] ~ dnorm(0, tau.eps[y])

            # Beta update
            # w[s, y] ~ dnorm(0, tau.w[s])
            beta[s, y+1] = beta[s, y] + zeta[s, y]

            # Lambda update
            lambda[s, y+1] = lambda[s, y] + eps[s, y] # winter to summer ratio

            # State update
            stt_s[s, y+1] = stt_s[s, y] + beta[s, y+1]# + w[s, y]
            stt_w[s, y+1] = stt_s[s, y+1] + lambda[s, y+1]

        }

    }

    # model for summer/winter
    for(i in 1:N){
        summer[i] ~ dbern(p[i])
        logit(p[i]) = inprod(B[site[i],], X[i,])
        winter[i] = 1 - summer[i]
    }

    ### likelihood ###

    for(i in 1:N){

        # State process
        mu_t[i] = stt_s[site[i], year[i]]*summer[i] + stt_w[site[i], year[i]]*winter[i]

        # observation process
        obs[i] ~ dnorm(mu_t[i], tau.alpha[site[i]]*summer[i] + tau.e[site[i]]*winter[i])

        count[i] ~ dpois(exp(obs[i]))

    }

}
