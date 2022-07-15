# This model considers two potentially unobserved states: summer and winter

model {

    ### PRIORS: ###

    # tau.w.hyp ~ dscaled.gamma(10, 2)

    for(s in 1:nsites){

        # Priors for initial states
        stt_s[s, 1] ~ dnorm(0, 1/10)       # prior for initial log population size
        beta[s, 1] ~ dnorm(0, 1/10)        # prior for initial log growth rate
        lambda[s, 1] ~ dnorm(0, 1/10)      # prior for log summer to winter ratio

        # # Prior for variation in abundance yearly slow change
        # tau.zeta[s] ~ dscaled.gamma(10, 2)    # same as Student-t with 2 degrees of freedom and standard deviation 10 for sig.zeta!
        # sig.zeta[s] = sqrt(1/tau.zeta[s])

        # # Prior for variation in winter ratio yearly slow change
        # tau.eps[s] ~ dscaled.gamma(10, 2)     # same as Student-t with 2 degrees of freedom and standard deviation 10 for sig.eps!
        # sig.eps[s] = sqrt(1/tau.eps[s])

        # Priors for observation error
        # tau.alpha[s] ~ dscaled.gamma(10, 10)   # summer
        # tau.e[s] ~ dscaled.gamma(10, 10)       # winter

        tau.alpha[s] ~ dgamma(2, 2)   # summer
        tau.e[s] ~ dgamma(2, 2)       # winter

        sig.alpha[s] = sqrt(1/tau.alpha[s])
        sig.e[s] = sqrt(1/tau.e[s])

        # Priors for random yearly abundance variation
        # tau.w[s] ~ dscaled.gamma(tau.w.hyp, 2)
        # sig.w[s] = sqrt(1/tau.w[s])

        # Prior for season covariate coefficients
        for(k in 1:K){
            B[s, k] ~ dnorm(0, 1)
        }

    }

    # Define variance-covariance matrices potentially based on distance between sites
    # for(i in 1:nsites){
    #
    #     mu.zeta[i] = 0
    #     mu.eps[i] = 0
    #
    #     for(j in 1:nsites){
    #         Omega.zeta[i, j] = tau.zeta[i]*tau.zeta[j]*equals(i, j)
    #         Omega.eps[i, j] = tau.eps[i]*tau.eps[j]*equals(i, j)
    #     }
    # }

    # likelihood

    # State updates
    for(s in 1:nsites){
        stt_w[s, 1] = stt_s[s, 1] + lambda[s, 1]
    }

    for (y in 1:(nyears-1)){

        # # Sample zeta and eps
        # zeta[1:nsites, y] ~ dmnorm(mu.zeta[1:nsites], Omega.zeta[1:nsites, 1:nsites])
        # eps[1:nsites, y] ~ dmnorm(mu.eps[1:nsites], Omega.eps[1:nsites, 1:nsites])

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
            stt_s[s, y+1] = stt_s[s, y] + beta[s, y]# + w[s, y]
            stt_w[s, y+1] = stt_s[s, y+1] + lambda[s, y+1]

        }

    }


    # model for summer/winter
    for(i in 1:N){
        summer[i] ~ dbern(p[i])
        logit(p[i]) = inprod(B[site[i],], X[i,])
        winter[i] = 1 - summer[i]
    }

    for(i in 1:N){

        # State process
        mu_t[i] = stt_s[site[i], year[i]]*summer[i] + stt_w[site[i], year[i]]*winter[i]

        # observation process
        obs[i] ~ dnorm(mu_t[i], tau.alpha[site[i]]*summer[i] + tau.e[site[i]]*winter[i])

        count[i] ~ dpois(exp(obs[i]))

    }

}
