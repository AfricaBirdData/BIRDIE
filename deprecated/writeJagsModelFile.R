#' Write a JAGS model to a temporary text file
#'
#' @return This function dumps the model below to a temporary file that calls
#' to rjags can use. It is meant to be a default model that is installed with
#' the package.
#' @export
#'
#' @examples
#' writeJagsModelFile()
writeJagsModelFile <- function(){

    filename <- tempfile(fileext = ".jags")

   mod <-  cat("

    # This model considers two time series: one for summer and one for winter
    # And a dynamic rate of chage (trend)

    model {

        # Priors
        mu_t[1] ~ dunif(log(1), log(1500)) # prior for initial population size
        beta[1] ~ dnorm(log(1), log(2)) # prior for initial long-term growth rate
        lambda[1] ~ dnorm(0, 0.001)

        tau.zeta ~ dunif(log(1), 1/log(3))
        tau.w ~ dunif(log(1), 1/log(10))
        tau.eps ~ dunif(log(1), 1/log(3))
        tau.alpha ~ dunif(log(1), 1/log(10))
        tau.e ~ dunif(log(1), 1/log(10))

        sig.zeta = 1/tau.zeta
        sig.w = 1/tau.w
        sig.eps = 1/tau.eps
        sig.alpha = 1/tau.alpha
        sig.e = 1/tau.e


        # likelihood

        # state process

        # summer
        for (t in 1:(N-1)){
            w[t] ~ dnorm(0, tau.w)
            zeta[t] ~ dnorm(0, tau.zeta)
            beta[t+1] = beta[t] + zeta[t]
            mu_t[t+1] = mu_t[t] + beta[t] + w[t]
        }

        # winter
        mu_wt[1] = mu_t[1] + lambda[1] # Initial winter state

        for (t in 1:(N-1)){
            eps[t] ~ dnorm(0, tau.eps)
            lambda[t+1] = lambda[t] + eps[t] # winter to summer ratio
            mu_wt[t+1] = mu_t[t+1] + lambda[t+1]
        }

        # observation process

        # summer
        for(s in 1:N){
            summer[s] ~ dnorm(mu_t[s], tau.alpha) # Summer counts
        }

        # winter
        for(w in 1:N){
            winter[w] ~ dnorm(mu_wt[w], tau.e) # Winter counts
        }

    }

",
        file = filename)

   return(filename)
}
