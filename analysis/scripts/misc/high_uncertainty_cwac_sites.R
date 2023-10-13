# This is a development version of a second round of selecting species to show on website

library(dplyr)
library(ggplot2)
library(BIRDIE)


rm(list = ls())


# Read in data
abu_data <- read.csv("analysis/out_nosync/ssm_pred_93_22_ZA_all_all.csv")

abu_data <- abu_data %>%
    mutate(excess = (summer.ci.upper - summer.est + 1)/(summer.est + 1))

hist(abu_data$excess)



# Calculate excess as max upper - max count ---------------------------------

site_data <- abu_data %>%
    filter(!is.na(summer.est)) %>%
    mutate(site_sp = paste0(site, "_", species)) %>%
    group_by(site_sp) %>%
    summarise(max_count = max(summer.count, na.rm = TRUE),
              max_up = max(summer.ci.upper, na.rm = TRUE)) %>%
    mutate(excess = (max_up - max_count)/max_count )

x <- c(seq(0, 0.9, length.out = 10), 0.95, 0.975, 0.99)
quantile(site_data$excess, probs = x, na.rm = TRUE)

plot(x, quantile(site_data$excess, probs = x, na.rm = TRUE))

high_excess <- site_data %>%
    filter(excess > quantile(excess, probs = 0.95, na.rm = TRUE))

low_excess <- site_data %>%
    filter(excess <= quantile(excess, probs = 0.95, na.rm = TRUE))

na_excess <- site_data %>%
    filter(is.na(excess))

nrow(high_excess) + nrow(low_excess) + nrow(na_excess)

hist(low_excess$excess)

# Identify site-species pairs in high excess
he_site_sp <- high_excess %>%
    distinct(site_sp)

# Number of sites in data with estimates
est_sites_data <- abu_data %>%
    filter(!is.na(excess))

length(unique(est_sites_data$site))


# Low excess data with estimates
le_sites_data <- abu_data %>%
    filter(!is.na(excess)) %>%
    mutate(site_sp = paste0(site, "_", species)) %>%
    filter(!site_sp %in% unique(he_site_sp$site_sp))

# Number of sites in low excess data with estimates
length(unique(le_sites_data$site))



# High excess data with estimates
he_sites_data <- abu_data %>%
    filter(!is.na(excess)) %>%
    mutate(site_sp = paste0(site, "_", species)) %>%
    filter(site_sp %in% unique(he_site_sp$site_sp))

he_sites_data %>%
    arrange(desc(excess)) %>%
    filter(site_sp == unique(site_sp)[1]) %>% {
        ggplot(.) +
            geom_point(aes(x = year, y = summer.count)) +
            geom_line(aes(x = year, y = summer.est)) +
            geom_line(aes(x = year, y = summer.ci.upper), col = "red") +
            ggtitle(paste(.$site_sp[1]))
    }

he_sites_data %>%
    arrange(excess) %>%
    filter(site_sp == unique(site_sp)[1]) %>% {
        ggplot(.) +
            geom_point(aes(x = year, y = summer.count)) +
            geom_line(aes(x = year, y = summer.est)) +
            geom_line(aes(x = year, y = summer.ci.upper), col = "red") +
            ggtitle(paste(.$site_sp[1]))
    }



le_sites_data %>%
    arrange(desc(excess)) %>%
    filter(site_sp == unique(site_sp)[1]) %>% {
        ggplot(.) +
            geom_point(aes(x = year, y = summer.count)) +
            geom_line(aes(x = year, y = summer.est)) +
            geom_line(aes(x = year, y = summer.ci.upper), col = "red") +
            ggtitle(paste(.$site_sp[1]))
    }



# Filter according to excess ----------------------------------------------

x <- c(seq(0, 0.9, length.out = 10), 0.95, 0.975, 0.99)
quantile(abu_data$excess, probs = x, na.rm = TRUE)

plot(x, quantile(abu_data$excess, probs = x, na.rm = TRUE))

high_excess <- abu_data %>%
    filter(excess > quantile(excess, probs = 0.95, na.rm = TRUE))

low_excess <- abu_data %>%
    filter(excess <= quantile(excess, probs = 0.95, na.rm = TRUE))

na_excess <- abu_data %>%
    filter(is.na(excess))

nrow(high_excess) + nrow(low_excess) + nrow(na_excess)

hist(low_excess$excess)

# Identify site-species pairs in high excess
he_site_sp <- high_excess %>%
    distinct(site, species) %>%
    mutate(site_sp = paste0(site, "_", species))

# Number of sites in data with estimates
est_sites_data <- abu_data %>%
    filter(!is.na(excess))

length(unique(est_sites_data$site))


# Low excess data with estimates
le_sites_data <- abu_data %>%
    filter(!is.na(excess)) %>%
    mutate(site_sp = paste0(site, "_", species)) %>%
    filter(!site_sp %in% unique(he_site_sp$site_sp))

# Number of sites in low excess data with estimates
length(unique(le_sites_data$site))



# High excess data with estimates
he_sites_data <- abu_data %>%
    filter(!is.na(excess)) %>%
    mutate(site_sp = paste0(site, "_", species)) %>%
    filter(site_sp %in% unique(he_site_sp$site_sp))

he_sites_data %>%
    arrange(desc(excess)) %>%
    filter(site_sp == unique(site_sp)[1]) %>% {
        ggplot(.) +
            geom_point(aes(x = year, y = summer.count)) +
            geom_line(aes(x = year, y = summer.est)) +
            geom_line(aes(x = year, y = summer.ci.upper), col = "red") +
            ggtitle(paste(.$site_sp[1]))
    }

he_sites_data %>%
    arrange(excess) %>%
    filter(site_sp == unique(site_sp)[1]) %>% {
        ggplot(.) +
            geom_point(aes(x = year, y = summer.count)) +
            geom_line(aes(x = year, y = summer.est)) +
            geom_line(aes(x = year, y = summer.ci.upper), col = "red") +
            ggtitle(paste(.$site_sp[1]))
    }



le_sites_data %>%
    arrange(desc(excess)) %>%
    filter(site_sp == unique(site_sp)[1]) %>% {
        ggplot(.) +
            geom_point(aes(x = year, y = summer.count)) +
            geom_line(aes(x = year, y = summer.est)) +
            geom_line(aes(x = year, y = summer.ci.upper), col = "red") +
            ggtitle(paste(.$site_sp[1]))
    }



# Filter according to last year counted -----------------------------------

last_year <- 2012

abu_data %>%
    group_by(year) %>%
    summarise(mean_excess = mean(excess, na.rm = TRUE)) %>%
    ungroup() %>%
    ggplot() +
    geom_col(aes(x = year, y = mean_excess))

high_excess <- abu_data %>%
    filter(!is.na(summer.count)) %>%
    group_by(site, species) %>%
    mutate(last_survey = max(year)) %>%
    ungroup() %>%
    filter(last_survey < last_year)

low_excess <- abu_data %>%
    filter(!is.na(summer.count)) %>%
    group_by(site, species) %>%
    mutate(last_survey = max(year)) %>%
    ungroup() %>%
    filter(last_survey >= last_year)

na_excess <- abu_data %>%
    filter(is.na(summer.count))

nrow(high_excess) + nrow(low_excess) + nrow(na_excess)

hist(low_excess$excess)

# Identify site-species pairs in high excess
he_site_sp <- high_excess %>%
    distinct(site, species) %>%
    mutate(site_sp = paste0(site, "_", species))

# Number of sites in data with estimates
est_sites_data <- abu_data %>%
    filter(!is.na(excess))

length(unique(est_sites_data$site))


# Low excess data with estimates
le_sites_data <- abu_data %>%
    filter(!is.na(excess)) %>%
    mutate(site_sp = paste0(site, "_", species)) %>%
    filter(!site_sp %in% unique(he_site_sp$site_sp))

# Number of sites in low excess data with estimates
length(unique(le_sites_data$site))



# High excess data with estimates
he_sites_data <- abu_data %>%
    filter(!is.na(excess)) %>%
    mutate(site_sp = paste0(site, "_", species)) %>%
    filter(site_sp %in% unique(he_site_sp$site_sp))

he_sites_data %>%
    arrange(desc(excess)) %>%
    filter(site_sp == unique(site_sp)[1]) %>% {
        ggplot(.) +
            geom_point(aes(x = year, y = summer.count)) +
            geom_line(aes(x = year, y = summer.est)) +
            geom_line(aes(x = year, y = summer.ci.upper), col = "red") +
            ggtitle(paste(.$site_sp[1]))
    }


le_sites_data %>%
    arrange(desc(excess)) %>%
    filter(site_sp == unique(site_sp)[3]) %>% {
        ggplot(.) +
            geom_point(aes(x = year, y = summer.count)) +
            geom_line(aes(x = year, y = summer.est)) +
            geom_line(aes(x = year, y = summer.ci.upper), col = "red") +
            ggtitle(paste(.$site_sp[1]))
    }
