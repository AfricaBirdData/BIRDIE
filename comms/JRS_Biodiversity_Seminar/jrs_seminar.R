# 30-06-2021

# In this script we run some simulations to produce plots for the JRS seminar

library(tidyverse)
library(gganimate)

rm(list = ls())

set.seed(43562345)

# set.seed(452345)

# Create temporal frame for abundance -------------------------------------

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

abund_plot <- df %>%
    ggplot() +
    geom_line(aes(x = t, y = logn, col = "blue")) +
    geom_segment(aes(x = t, y = logn,
                     xend = t, yend = obs, group = seq_along(t),
                     col = "darkgrey"), linetype = "dashed", size = 1) +
    geom_point(aes(x = t, y = obs, group = seq_along(t)), size = 3) +
    geom_line(data = df[!is.na(df$obs),],
              aes(x = t, y = obs, col = "black")) +
    geom_point(data = df[!is.na(df$obs),],
               aes(x = t, y = logn, group = seq_along(t)),
               col = "darkgrey", shape = 1, size = 3, stroke = 1) +
    scale_colour_manual(values = c("black", "blue", "darkgrey"),
                        labels = c("Observed", "True", "Error"),
                        name = "") +
    xlab("Time") + ylab("Abundance") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          legend.text = element_text(size = 13),
          legend.position = "top")

abund_plot +
    transition_reveal(t)

ggsave(plot = abund_plot, filename = "comms/output/abund_obs.png")
anim_save("comms/output/abund_obs_anim.gif")


# Create temporal frame for abundance -------------------------------------

# Number of observations
N <- 100

# Preliminary time increments(months)
dt <- ceiling(rexp(N-1, 0.5))

# Times
t <- c(0, cumsum(dt))

# Create regular observation times
t_obs <- seq(20, max(t) - 20, length.out = 10)

# Integrate with the process simulation times
t <- sort(c(t[!t %in% t_obs], t_obs))

# Recalculate dt
dt <- lead(t) - t


# Create a time series of an occupancy process ----------------------------

set.seed(4526)

# Occupancy probabilities vector in logit
logitp <- vector("numeric", length = length(t))

# Initial probability (logit)
logitp[1] <- -1

# Diffusion standard deviation
sigma <- 0.01

# Change standard deviation
sig_beta <- 0.005

# Initial probability change
beta <- -0.01

# Generate abundance time series
for(i in seq_along(logitp[-1])){
    logitp[i+1] <- logitp[i] + beta + rnorm(1, 0, sigma)
    beta <- beta + rnorm(1, 0, sig_beta*sqrt(dt[i]))
}

plot(t, logitp, type = "l")
points(t_obs, logitp[t %in% t_obs])

# Calculate p and draw a sample at the time of observation
p <- exp(logitp) / (1 + exp(logitp))

occu <- rbinom(length(t_obs), 1, p[t %in% t_obs])

plot(t, p, type = "l", ylim = c(0, 1))
points(t_obs, occu)


# Create a time series of detections --------------------------------

# Probability of observing given presence
pcond <- 0.5

occu_obs <- rbinom(length(occu), 1, ifelse(occu == 1, 0.5, 0))

plot(t, p, type = "l", ylim = c(0, 1))
points(t_obs, occu)
points(t_obs, occu_obs, col = "red")

# Create animated plot ----------------------------------------------------

pres <- rep(NA, length = length(p))
obs <- rep(NA, length = length(p))
pres[t %in% t_obs] <- occu
obs[t %in% t_obs] <- occu_obs

# jitter obs a little
obs <- ifelse(obs == 1, 0.95, 0.05)

df <- data.frame(t = t, dt = dt, p = p, pres = pres, obs = obs)

occu_plot <- df %>%
    ggplot() +
    geom_line(aes(x = t, y = p), col = "blue") +
    geom_point(aes(x = t, y = pres, group = seq_along(t), col = "blue"), size = 3) +
    geom_point(aes(x = t, y = obs, group = seq_along(t), col = "black"), size = 3) +
    scale_colour_manual(values = c("black", "blue"),
                        labels = c("Observed", "True"),
                        name = "") +
    xlab("Time") + ylab("Probability of occurrence") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          legend.text = element_text(size = 13),
          legend.position = "top")

occu_plot +
    transition_reveal(t)

ggsave(plot = occu_plot, filename = "comms/output/occu_obs.png")
anim_save("comms/output/occu_obs_anim.gif")

