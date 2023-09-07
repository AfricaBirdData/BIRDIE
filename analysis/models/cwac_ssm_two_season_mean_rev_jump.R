# This model considers two states: abundance in summer and winter
# It estimates these states for a number of sites, using counts taken at those sites
# Some of the sites may be missing counts for some of the years.

# This is a mean reverting process in which the mean is estimated as a function
# of covariates

# NOTE THAT THE MODEL ASSUMES COUNTS ARE ADDED ONE TO AVOID ZERO COUNTS AND NON-IDENTIFIABILITY OF PARAMETERS.
# SEE RE-DEFINITION OF MU_T_1 INTO MU_T AT THE END OF THE MODEL

model {


    # model for summer/winter
    for(i in 1:N){
        winter[i] = 1 - summer[i]
    }

    ### PRIORS: ###

    for(s in 1:nsites){

        # Priors for initial states
        mu_ini[s] ~ dnorm(psi_s_sc[s, 1], 1)              # prior for initial log abundance
        beta_ini[s] ~ dnorm(0, 1)
        lambda_ini[s] ~ dnorm(psi_w_sc[s, 1] - mu_ini[s], 1)             # prior for initial log summer to winter ratio

        phi[s] ~ dbeta(2, 5)
        #phi.lambda[s] ~ dbeta(2, 2)

        # tau.alpha[s] ~ dscaled.gamma(1, 4)   # same as Student-t with sd standard deviation and df degrees of freedom where dscaled.gamma(sd, df) for sigma
        # tau.e[s] ~ dscaled.gamma(1, 4)       # same as Student-t with sd standard deviation and df degrees of freedom where dscaled.gamma(sd, df) for sigma

        tau.alpha[s] ~ dgamma(1, 1)   # Using gamma to prevent sigma to get too close to zero
        tau.e[s] ~ dgamma(1, 1)       # Using gamma to prevent sigma to get too close to zero

        sig.alpha[s] = sqrt(1/tau.alpha[s])   # observer error standard deviation summer
        sig.e[s] = sqrt(1/tau.e[s])           # observer error standard deviation winter

        psi_s_sc0[s] ~ dnorm(med_log_count[s], 1)                # Prior for site-level intercept for summer abundance
        psi_w_sc0[s] ~ dnorm(med_log_count[s] + med_ws[s], 1)               # Prior for site-level intercept for winter abundance

        psi_s_sc[s, 1] = psi_s_sc0[s]
        psi_w_sc[s, 1] = psi_w_sc0[s]

        pj[s] ~ dbeta(1, 50)

    }

    tau.jump ~ dgamma(5, 2)

    # These covariates are used to model the mean abundance over time.
    # for(m in 1:M){
    #
    #     mu_g[m] ~ dnorm(0, 1)
    #     sig_g[m] ~ dnorm(1, 5) T(0,)
    #
    #     for(s in 1:nsites){
    #         G[s, m] ~ dnorm(mu_g[m], sig_g[m])
    #     }
    #
    # }

    # Priors for stochastic summer population changes
    for(s in 1:nsites){
        tau.zeta[s] ~ dscaled.gamma(1, 20)   # same as Student-t with sd standard deviation and df degrees of freedom where dscaled.gamma(sd, df) for sigma
        # tau.zeta[y] ~ dgamma(6, 3)
        sig.zeta[s] = sqrt(1/tau.zeta[s])
        # zeta[s, 1] ~ dnorm(0, tau.zeta[s])
        for(y in 1:nyears){
            Ij_s[s, y] ~ dbern(pj[s])
            z_s[s, y] ~ dnorm(0, tau.jump)
            zeta[s, y] ~ dnorm(0, tau.zeta[s])
        }
    }

    # Priors for changes in winter-to-summer abundance ratio
    for(s in 1:nsites){
        tau.eps[s] ~ dscaled.gamma(1, 20)   # same as Student-t with sd standard deviation and df degrees of freedom where dscaled.gamma(sd, df) for sigma
        # tau.eps[y] ~ dgamma(6, 3)
        sig.eps[s] = sqrt(1/tau.eps[s])
        for(y in 1:nyears){
            Ij_w[s, y] ~ dbern(pj[s])
            z_w[s, y] ~ dnorm(0, tau.jump)
            eps[s, y] ~ dnorm(0, tau.eps[s])
        }
    }

    ### LATENT STATES ###

    # Estimate covariate effects on abundance
    # for(y in 1:nyears){
    #     for(s in 1:nsites){
    #         eta[s, y] = inprod(G[s, ], U[y,,s])
    #     }
    # }

    for(s in 1:nsites){

        # Initial Change in summer abundance
        beta_p[s, 1] = beta_ini[s]    # this is used later, and it is necessary to model a smooth beta trajectory
        jump_s[s, 1] = Ij_s[s, 1]*z_s[s, 1]
        stt_s_sc[s, 1] = mu_ini[s] + beta_ini[s] + jump_s[s, 1]

        # Initial Winter-to-summer abundance ratios
        lambda_p[s, 1] = lambda_ini[s]    # this is used later, and it is necessary to model a smooth beta trajectory
        jump_w[s, 1] = Ij_w[s, 1]*z_w[s, 1]
        stt_w_sc[s, 1] = stt_s_sc[s, 1] + lambda_ini[s]
    }

    for(s in 1:nsites){
        for(y in 1:(nyears-1)){

            # Change in summer abundance
            zeta_p[s, y] = psi_s_sc[s, y] - stt_s_sc[s, y]  # Expected mean-reverting change in summer
            beta_p[s, y+1] = beta_p[s, y] + phi[s]*(zeta_p[s, y] - beta_p[s, y]) + zeta[s, y]
            jump_s[s, y+1] = Ij_s[s, y+1]*z_s[s, y+1]

            # Summer state update
            beta[s, y] = beta_p[s, y+1] + jump_s[s, y+1]
            psi_s_sc[s, y+1] = psi_s_sc[s, y] + jump_s[s, y+1]
            stt_s_sc[s, y+1] = stt_s_sc[s, y] + beta[s, y]

            # Winter-to-summer abundance ratios
            eps_p[s, y] = psi_w_sc[s, y] - stt_s_sc[s, y+1] # Expected mean-reverting change in winter
            lambda_p[s, y+1] = lambda_p[s, y] + phi[s]*(eps_p[s, y] - lambda_p[s, y]) + eps[s, y+1]
            jump_w[s, y+1] = Ij_w[s, y+1]*z_w[s, y+1]

            # Winter state update
            lambda[s, y] = lambda_p[s, y+1] + jump_w[s, y+1]
            psi_w_sc[s, y+1] = psi_w_sc[s, y] + jump_s[s, y+1]
            stt_w_sc[s, y+1] = stt_s_sc[s, y+1] + lambda[s, y]

        }
    }

    # Re-scale to original site scale
    for(s in 1:nsites){
        for(y in 1:nyears){
            psi_s[s, y] = psi_s_sc[s, y] + log(sc_site[s])
            psi_w[s, y] = psi_w_sc[s, y] + log(sc_site[s])
            stt_s[s, y] = stt_s_sc[s, y] + log(sc_site[s])
            stt_w[s, y] = stt_w_sc[s, y] + log(sc_site[s])
        }
    }

    # Latent process
    for(i in 1:N){

        # State process
        mu_t[i] = stt_s_sc[site[i], year[i]]*summer[i] + stt_w_sc[site[i], year[i]]*winter[i]
        tau_t[i] = tau.alpha[site[i]]*summer[i] + tau.e[site[i]]*winter[i]

        # observation process. Should this be dt?
        count[i] ~ dnorm(mu_t[i], tau_t[i])

        #count[i] ~ dpois(exp(obs[i]))

    }


    ### LIKELIHOOD ###

    for(i in 1:N_data_obs){

        # Log likelihood
        # log.lik[i] = logdensity.pois(count[obs_idx[i]], exp(obs[obs_idx[i]]))
        log.lik[i] = logdensity.norm(count[obs_idx[i]], mu_t[obs_idx[i]], tau_t[obs_idx[i]])

    }

    # Prepare mu_t_1 for export: SUBTRACT 1 BECAUSE WE ADDED 1 TO ALL COUNTS
    # BETAS AND LAMBDAS SHOULD BE UNAFFECTED BY THIS
    # for(i in 1:N){
    #     mu_t[i] = log(max(0.001, exp(mu_t_1[i]) - 1)*sc_site[site[i]])
    # }


}
