library(BIRDIE)
library(dplyr)

rm(list = ls())

year <- 2006

config <- configPreambJAGS(year = year, server = TRUE)
species <- c(4, 6)

# cwacsites <- CWAC::listCwacSites(.region_type = "province", .region = "northern cape")
#
# cwacsites %>%
#     filter(LocationName == "Du Toit's Pan")

sites <- c(26352535, 28462448)

for(s in seq_along(sites)){

    site <- sites[s]

    for(i in seq_along(species)){

        config <- configPreambJAGS(year = year, server = FALSE)

        sp_code <- species[i]

        # Species name
        sp_name <- BIRDIE::barberspan %>%
            dplyr::filter(SppRef == sp_code) %>%
            mutate(name = paste(Common_species, Common_group)) %>%
            pull(name) %>%
            unique()

        print(paste0("Working on species ", sp_code, " (", i, " of ", length(config$species), ")"))

        ppl_run_pipe_abu(sp_code = sp_code,
                         site = site,
                         year = year,
                         config = config,
                         steps = c("fit", "summ"))

    }

    # Estimate species type indicators
    ppl_run_pipe_abu(sp_code = species,
                     site = site,
                     year = year,
                     config = config,
                     steps = c("fit", "summ"))
}

library(dplyr)
library(ggplot2)

sp1 <- read.csv("comms/workshop_feb_22/ssm_pred_26352535_93_19_4.csv")
sp2 <- read.csv("comms/workshop_feb_22/ssm_pred_26352535_93_19_6.csv")
group <- read.csv("comms/workshop_feb_22/ssm_pred_26352535_93_19_group.csv")

colors <- c("#71BD5E", "#B590C7")

ggsave(
    sp1 %>%
        ggplot() +
        geom_line(aes(x = year, y = summer.logest), linetype = 1) +
        geom_line(aes(x = year, y = summer.logest.ci.lower), linetype = 2) +
        geom_line(aes(x = year, y = summer.logest.ci.upper), linetype = 2) +
        geom_point(aes(x = year, y = log.summer.count)) +
        xlab("Year") + ylab("Abundance") +
        ggtitle("Great Crested Grebe") +
        theme_classic() +
        theme(text = element_text(size = 16))
    , filename = "comms/workshop_feb_22/plot_sp1.png")


ggsave(
    sp2 %>%
        ggplot() +
        geom_line(aes(x = year, y = summer.logest), linetype = 1) +
        geom_line(aes(x = year, y = summer.logest.ci.lower), linetype = 2) +
        geom_line(aes(x = year, y = summer.logest.ci.upper), linetype = 2) +
        geom_point(aes(x = year, y = log.summer.count)) +
        xlab("Year") + ylab("Abundance") +
        ggtitle("Little Grebe") +
        theme_classic() +
        theme(text = element_text(size = 16))
    , filename = "comms/workshop_feb_22/plot_sp2.png")


ggsave(
    group %>%
        ggplot() +
        geom_line(aes(x = year, y = summer.logest), linetype = 1) +
        geom_line(aes(x = year, y = summer.logest.ci.lower), linetype = 2) +
        geom_line(aes(x = year, y = summer.logest.ci.upper), linetype = 2) +
        geom_point(aes(x = year, y = log.summer.count)) +
        xlab("Year") + ylab("Abundance") +
        ggtitle("Multispecies") +
        theme_classic() +
        theme(text = element_text(size = 16))
    , filename = "comms/workshop_feb_22/plot_group.png")


