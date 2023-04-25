library(dplyr)
library(ABAP)
library(CWAC)
library(ggplot2)

rm(list = ls())



# Summary of efforts in ABAP and CWAC -------------------------------------

# Download ABAP data for any species for the last ten years
abap_data <- getAbapData(.spp_code = 6,
                         .region_type = "country",
                         .region = "South Africa",
                         .years = 2011:2021)

# Download South African ABAP pentads
pentads_sa <- getRegionPentads(.region_type = "country",
                               .region = "South Africa")

# Calculate number of hours per pentad
spt_effort <- abap_data %>%
    group_by(Pentad) %>%
    summarise(hours = sum(TotalHours))

# Plot spatial effort
pentads_sa %>%
    left_join(spt_effort, by = c("pentad" = "Pentad")) %>%
    ggplot() +
    geom_sf(aes(fill = hours/1000), lwd = NA) +
    scale_fill_viridis_c(name = "khours", option = "B")

pentads_sa %>%
    left_join(spt_effort, by = c("pentad" = "Pentad")) %>%
    ggplot() +
    geom_sf(aes(fill = log(hours+1)), lwd = NA) +
    scale_fill_viridis_c(name = "log(hours)", option = "B", direction = 1)

# Calculate effort over time
temp_effort <- abap_data %>%
    mutate(year = lubridate::year(StartDate)) %>%
    group_by(year) %>%
    summarise(hours = sum(TotalHours))

# Plot temporal efforts
temp_effort %>%
    ggplot(aes(x = year, y = hours/1000)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    ylab("Hours / 1000")


# Download CWAC sites
cwac_sites <- listCwacSites(.region_type = "country",
                            .region = "South Africa")

cwac_counts <- getCwacSppCounts(6)

cwac_sites <- cwac_sites %>%
    mutate(lastyear = lubridate::year(LastSurvey),
           firstyear = lubridate::year(FirstSurvey))

# Remove NA first and last year
cwac_sites <- cwac_sites %>%
    filter(!is.na(firstyear),
           !is.na(lastyear),
           firstyear > 1990)

# Filter South African sites
cwac_sites <- cwac_sites %>%
    filter(Y < -10)

# Calculate duration and filter out variables
cwac_sites <- cwac_sites %>%
    mutate(duration = lastyear - firstyear) %>%
    arrange(firstyear, duration) %>%
    select(LocationCode, X, Y, firstyear, lastyear, duration) %>%
    distinct() %>%
    mutate(LocationCode = factor(LocationCode, levels = .$LocationCode))

# Plot timeline
cwac_sites %>%
    ggplot() +
    geom_segment(aes(x = firstyear, xend = lastyear,
                     y = LocationCode, yend = LocationCode, col = duration), size = 0.5) +
    scale_color_viridis_c(name = "dur (years)") +
    xlab("Year") + ylab("Site") +
    ggtitle("Timeline of site counts (each fine horizontal line is a site)") +
    theme(axis.text.y = element_blank())

# Plot start of monitoring
cwac_sites %>%
    count(firstyear) %>%
    ggplot() +
    geom_point(aes(x = firstyear, y = n))

# Plot end of monitoring
cwac_sites %>%
    filter(lastyear < 2022) %>%
    count(lastyear) %>%
    ggplot() +
    geom_point(aes(x = lastyear, y = n))

# Plot duration
cwac_sites %>%
    ggplot() +
    geom_histogram(aes(x = duration))




# CWAC duration mixture test ----------------------------------------------

alpha <- c(1, 10, 25)
beta <- c(1, 1, 1)
mixmod <- mixtools::gammamixEM(x = cwac_sites$duration + 1, alpha = alpha, beta = beta,
                               k = 3, maxit = 2000)
mixmod$loglik

plot(mixmod$all.loglik)
mixmod$gamma.pars
mixmod$lambda

mu_fit <- mixmod$gamma.pars[1,] * mixmod$gamma.pars[2,]

hist(cwac_sites$duration + 1, freq = F)
curve(dgamma(x, shape = mixmod$gamma.pars[1,1], scale = mixmod$gamma.pars[2,1])*mixmod$lambda[1], col = "red", add = T)
curve(dgamma(x, shape = mixmod$gamma.pars[1,2], scale = mixmod$gamma.pars[2,2])*mixmod$lambda[2], add = T)
curve(dgamma(x, shape = mixmod$gamma.pars[1,3], scale = mixmod$gamma.pars[2,3])*mixmod$lambda[3], add = T)
