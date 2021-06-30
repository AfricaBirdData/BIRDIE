# 30-06-2021

# In this script we run some simulations to produce plots for the JRS seminar

library(tidyverse)
library(gganimate)

rm(list = ls())

set.seed(43562345)

# set.seed(452345)

# Create temporal frame ---------------------------------------------------

# Number of observations
N <- 100

# Preliminary time increments(months)
dt <- ceiling(rexp(N-1, 0.5))

# Times
t <- c(0, cumsum(dt))

# Create regular observation times
t_obs <- seq(20, max(t) - 50, length.out = 5)

# Integrate with the process simulation times
t <- sort(c(t[!t %in% t_obs], t_obs))

# Recalculate dt
dt <- lead(t) - t


# Create a time series of an abundance process ----------------------------

# Abundance vector in log scale
logn <- vector("numeric", length = length(t))

# Initial population size (log)
logn[1] <- 5

# Diffusion standard deviation
sigma <- 0.1

# Change standard deviation
sig_beta <- 0.05

# Initial population change
beta <- 0

# Generate abundance time series
for(i in seq_along(logn[-1])){
    logn[i+1] <- logn[i] + beta + rnorm(1, 0, sigma)
    beta <- beta + rnorm(1, 0, sig_beta*sqrt(dt[i]))
}

plot(t, logn, type = "l")
points(t_obs, logn[t %in% t_obs])


# Create a time series of abundance observations --------------------------

# Observation error variance
sigma_obs <- 1.5

logn_obs <- logn[t %in% t_obs] + rnorm(length(t_obs), 0, sigma_obs)

plot(t, logn, type = "l")
points(t_obs, logn_obs)


# Create animated plot ----------------------------------------------------

obs <- rep(NA, length = length(logn))
obs[t %in% t_obs] <- logn_obs

df <- data.frame(t = t, dt = dt, logn = logn, obs = obs)

df %>%
    ggplot() +
    geom_line(aes(x = t, y = logn), col = "blue") +
    geom_segment(aes(x = t, y = logn, xend = t, yend = obs, group = seq_along(t)),
                 col = "darkgrey", linetype = "dashed", size = 1) +
    geom_point(aes(x = t, y = obs, group = seq_along(t)), size = 3) +
    geom_line(data = df[!is.na(df$obs),], aes(x = t, y = obs)) +
    xlab("Time") + ylab("Abundance") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 13)) +
    transition_reveal(t)

anim_save("comms/output/abund_obs_anim.gif")
