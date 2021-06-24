
library(BIRDIE)

rm(list = ls())


# DATA --------------------------------------------------------------------

counts <- barberspan
counts <- prepCtSsmData(counts, species = NULL)

obs = log(counts$count + 0.1)
summer = counts$summer
winter = counts$winter
dsummer = counts$dseason[counts$summer == 1]
dwinter = counts$dseason[counts$winter == 1]
tgt_idx = counts$tgt_idx
dt = counts$dt
N = nrow(counts)
Ns = sum(counts$summer)
Nw = sum(counts$winter)

# PRIORS -----------------------------------------------------------------

tau.zeta = 0.5
tau.eps = 0.5
sig.B = 0.1
tau.alpha = 0.1
tau.e = 0.1
tau.o = 0.1
gamma = 0.1

sig.zeta = 1/tau.zeta
sig.eps = 1/tau.eps
sig.alpha = 1/tau.alpha
sig.e = 1/tau.e
sig.o = 1/tau.o

beta <- vector("numeric", length = Ns)
zeta <- vector("numeric", length = Ns-1)

lambda <- vector("numeric", length = Nw)
eps <- vector("numeric", length = Nw-1)

theta <- vector("numeric", length = N)
mu <- obs


# MODEL FOR POPULATION INCREMENTS ----------------------------------------

# model for summer population increments
beta[1] = 0.1 # prior for initial long-term growth rate

for(i in 1:(Ns-1)){

    zeta[i] = rnorm(1, 0, sig.zeta*sqrt(dsummer[i]))
    beta[i+1] = beta[i] + zeta[i]

}

plot(beta)
plot(zeta)

# model for winter population increments
lambda[1] = 0.1

for(i in 1:(Nw-1)){

    eps[i] = rnorm(1, 0, sig.eps*sqrt(dwinter[i]))
    lambda[i+1] = lambda[i] + eps[i]

}

plot(lambda)
plot(eps)


# MODEL FOR TARGET LATENT STATES -----------------------------------------

# indices to keep track of the different seasons
s <- vector("integer", length = Ns)
w <- vector("integer", length = Nw)

s[1] = 1
w[1] = 1

theta[1] = mu[1] # this is actually irrelevant, but needs to be filled in

for(t in 1:(N-1)){

    theta[t+1] = theta[t] - theta[t] * (winter[t] + summer[t]) + # if no winter or summer
        (mu[t] + beta[s[t]]) * summer[t] + # if summer
        (mu[t] + lambda[w[t]]) * winter[t] # if winter

    s[t+1] = s[t] + 1 * summer[t] * (Ns > s[t]) # if summer
    w[t+1] = w[t] + 1 * winter[t] * (Nw > w[t]) # if winter

}


# MODEL FOR POPULATION LATENT STATE ---------------------------------------

mu[1] ~ dunif(log(1), log(1500)) # prior for initial population size

for(t in 1:(N-1)){

    Diff[t] = (sig.B/(2*gamma))*(1 - exp(-2*gamma*dt[t]))
    B[t] ~ dnorm(0, 1/sqrt(Diff[t]))
    mu[t+1] = mu[t]*exp(-gamma*dt[t]) + theta[tgt_idx[t]]*(1-exp(-gamma*dt[t])) + B[t]

}


# MODEL FOR OBSERVATION PROCESS -----------------------------------------

for(t in 1:N){
    obs[t] ~ dnorm(mu[t], tau.alpha * summer[t] +
                       tau.e * winter[t] +
                       tau.o - tau.o * (summer[t] + winter[t]))
}
