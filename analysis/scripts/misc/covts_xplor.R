library(BIRDIE)
library(tidyverse)

load("analysis/data/Barberspan_NEX-GDPP_clim_covts.rda")

dat$Date <- lubridate::date(dat$Date)

summary(dat)

dat %>%
    slice(1:4) %>%
    ggplot() +
    geom_point(aes(x = Date, y = Precip_mean))
