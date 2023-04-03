# This model considers two potentially unobserved states: abundance in summer and winter
# It estimates these states for a number of sites, using counts taken at those sites
# Some of the sites may be missing counts for some of the years.
# State changes are estimated in a hierarchical way such that changes across sites
# in a given year come from a common distribution.

# NOTE THAT THE MODEL ASSUMES COUNTS ARE ADDED ONE TO AVOID ZERO COUNTS AND NON-IDENTIFIABILITY OF PARAMETERS.
# SEE RE-DEFINITION OF MU_T_1 INTO MU_T AT THE END OF THE MODEL

model {

    ### PRIORS: ###

    for(s in 1:nsites){

        # Priors for initial states
        zeta_ini[s] ~ dnorm(0, 3)              # prior for initial extra log abundance uncertainty
        lambda_ini[s] ~ dnorm(0, 3)             # prior for log summer to winter ratio extra uncertainty
        stt_w[s, 1] = stt_s[s, 1] + lambda[s, 1]  # This is not a prior but inherits directly from two priors

        phi.mu[s] ~ dbeta(2, 2)
        phi.lambda[s] ~ dbeta(2, 2)

        tau.alpha[s] ~ dgamma(2, 0.5)
        tau.e[s] ~ dgamma(2, 0.5)

        sig.alpha[s] = sqrt(1/tau.alpha[s])   # observer error standard deviation summer
        sig.e[s] = sqrt(1/tau.e[s])           # observer error standard deviation winter

        # Prior for season covariate coefficients
        # for(k in 1:K){
        #     B[s, k] ~ dnorm(0, 1)
        # }

        # Prior for expected population change coefficients
        # G[s, 1] ~ dnorm(0, 0.5)

        for(m in 1:M){
            G[s, m] ~ dnorm(0, 1)
        }
    }

    # Priors for stochastic summer population changes
    for(y in 1:nyears){
        tau.zeta[y] ~ dscaled.gamma(3, 4)    # same as Student-t with 4 degrees of freedom and standard deviation 1 for sig.zeta!
        # tau.zeta[y] ~ dgamma(6, 3)
        sig.zeta[y] = sqrt(1/tau.zeta[y])
        for(s in 1:nsites){
            zeta[s, y] ~ dnorm(0, tau.zeta[y])
        }
    }

    # Priors for changes in winter-to-summer abundance ratio
    for(y in 1:nyears){
        tau.eps[y] ~ dscaled.gamma(3, 4)    # same as Student-t with 4 degrees of freedom and standard deviation 1 for sig.eps!
        # tau.eps[y] ~ dgamma(6, 3)
        sig.eps[y] = sqrt(1/tau.eps[y])
        for(s in 1:nsites){
            eps[s, y] ~ dnorm(0, tau.eps[y])
        }
    }

    ### LATENT STATES ###

    ## SUMMER ABUNDANCE

    # Estimate covariate effects on abundance
    for(y in 1:nyears){
        for(s in 1:nsites){
            eta[s, y] = inprod(G[s, ], U[y,,s])
        }
    }

    # Summer abundance estimates
    for(s in 1:nsites){
        mu_beta[s, 1] = eta[s, 1] + zeta[s, 1] + zeta_ini[s]
    }

    for(s in 1:nsites){
        for(y in 2:nyears){
            mu_beta[s, y] = eta[s, y] + zeta[s, y]
        }
    }

    ## WINTER ABUNDANCE

    # Note that we need more lambdas than betas. This
    # because we need to estimate winter to summer proportion for every year, but
    # we don't estimate population change for the last year.

    # Winter-to-summer abundance ratios
    for(s in 1:nsites){
        # Initial state
        lambda[s, 1] = eps[s, 1]  + lambda_ini[s]
    }

    for(y in 1:(nyears - 1)){
        for(s in 1:nsites){
            lambda[s, y+1] = phi.lambda[s]*lambda[s, y] + (1-phi.lambda[s])*eps[s, y]                          # winter to summer ratio
        }
    }

    # Update population states
    for(s in 1:nsites){
        # Initial state
        stt_s[s, 1] = mu_beta[s, 1]
    }

    for (y in 1:(nyears-1)){
        for(s in 1:nsites){

            # State update
            stt_s[s, y+1] = phi.mu[s]*stt_s[s, y] + (1-phi.mu[s])*mu_beta[s, y+1]
            stt_w[s, y+1] = stt_s[s, y+1] + lambda[s, y+1]

            beta[s, y] = stt_s[s, y+1] - stt_s[s, y]

        }
    }


    # model for summer/winter
    for(i in 1:N){
        # summer[i] ~ dbern(p[i])
        # logit(p[i]) = inprod(B[site[i],], X[i,])
        winter[i] = 1 - summer[i]
    }

    ### likelihood ###

    for(i in 1:N){

        # State process
        mu_t_1[i] = stt_s[site[i], year[i]]*summer[i] + stt_w[site[i], year[i]]*winter[i]

        # observation process
        obs[i] ~ dnorm(mu_t_1[i], tau.alpha[site[i]]*summer[i] + tau.e[site[i]]*winter[i])

        count[i] ~ dpois(exp(obs[i]))

    }

    # Prepare mu_t_1 for export: SUBTRACT 1 BECAUSE WE ADDED 1 TO ALL COUNTS
    # BETAS AND LAMBDAS SHOULD BE UNAFFECTED BY THIS
    for(i in 1:N){
        mu_t[i] = log(max(0.001, exp(mu_t_1[i]) - 1))
    }


}
