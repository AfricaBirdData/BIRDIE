library(BIRDIE)
library(tidyverse)

counts <- barberspan

# Prepare times and time increments
datetimes <- data.frame(startDate = unique(counts$startDate)) %>%
    dplyr::arrange(startDate) %>%
    dplyr::mutate(dt = as.numeric(difftime(dplyr::lead(startDate), startDate, units = "days")),
                  # dt = ifelse(is.na(dt), 0, dt),
                  t = c(0, cumsum(dt[!is.na(dt)])))

counts <- counts %>%
    dplyr::left_join(datetimes, by = "startDate")

counts %>%
    group_by(t) %>%
    summarize(n = sum(count),
              season = unique(Season)) %>%
    ggplot() +
    geom_point(aes(x = t, y = log(n), col = season)) +
    geom_path(aes(x = t, y = log(n)))

acf(log(counts$count))

counts %>%
    group_by(t) %>%
    summarize(n = sum(count),
              season = unique(Season)) %>%
    ggplot() +
    geom_point(aes(x = t, y = log(n))) +
    geom_path(aes(x = t, y = log(n))) +
    facet_wrap("season", nrow = 3)

counts %>%
    group_by(t) %>%
    summarize(date = unique(startDate),
              season = unique(Season),
              dt = unique(dt),
              n = sum(count)) %>%
    print(n = Inf)

cc <- counts %>%
    group_by(t) %>%
    summarize(date = unique(startDate),
              season = unique(Season),
              dt = unique(dt),
              n = sum(count)) %>%
    mutate(year = lubridate::year(date))

# Find index of target counts ("W" and "S")
cc <- cc %>%
    mutate(tgt_idx = findNextIndex(season, c("W", "S")))

# Create dseason, number of seasons between observations
dseason <- vector("integer", length = nrow(cc))
ss <- c(diff(cc$year[cc$season == "S"]), NA)
ww <- c(diff(cc$year[cc$season == "W"]), NA)
ids <- 1
idw <- 1

for(i in 1:nrow(cc)){

    # update dseason
    if(cc$season[i] == "S"){
        dseason[i] <- ss[ids]
        ids <- ids + 1
    }

    if(cc$season[i] == "W"){
        dseason[i] <- ww[idw]
        idw <- idw + 1
    }

}

cc <- cc %>%
    mutate(dseason = dseason)

print(cc, n = Inf)
