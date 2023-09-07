library(dplyr)
library(ggplot2)

dat_4 <- read.csv("analysis/out_nosync/4/abu_model_data_jagsUI_93_22_4.csv")

# Add one to all counts so that we can take logs
counts <- dat_4 %>%
    dplyr::mutate(count = count + 1)

# Scale response so that all sites are similar
sc_site <- counts %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(sd_count = stats::sd(count, na.rm = TRUE)) %>%
    dplyr::pull(sd_count)

counts <- counts %>%
    dplyr::group_by(site_id) %>%
    dplyr::mutate(sd_count = stats::sd(count, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(count = count / sd_count)

counts %>%
    ggplot(aes(x = year, y = count)) +
    geom_point() +
    geom_line() +
    facet_wrap("site_id")


dat_4 %>%
    ggplot() +
    geom_histogram(aes(x = count)) +
    facet_wrap("site_id")

dat_4 %>%
    ggplot() +
    geom_histogram(aes(x = count)) +
    facet_wrap("site_id")

dat_87 <- read.csv("analysis/out_nosync/87/abu_model_data_jagsUI_93_22_87.csv")

# Add one to all counts so that we can take logs
counts <- dat_87 %>%
    dplyr::mutate(count = count + 1)

# Scale response so that all sites are similar
sc_site <- counts %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(sd_count = stats::sd(count, na.rm = TRUE)) %>%
    dplyr::pull(sd_count)

counts <- counts %>%
    dplyr::group_by(site_id) %>%
    dplyr::mutate(sd_count = stats::sd(count, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(count = count / sd_count)

counts %>%
    ggplot(aes(x = year, y = count)) +
    geom_point() +
    geom_line() +
    facet_wrap("site_id")

dat_87 %>%
    ggplot() +
    geom_histogram(aes(x = count)) +
    facet_wrap("site_id")
