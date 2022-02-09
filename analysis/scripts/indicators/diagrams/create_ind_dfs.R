library(dplyr)
library(DiagrammeR)
library(DiagrammeRsvg)

rm(list = ls())


# Create element list -----------------------------------------------------

elem_list <- list(taxon = c("species", "sp_type"),
                  geo = c("site", "site_type", "global"),
                  time = "time",
                  cond = c("at", "across"),
                  indicator = c("distr", "abund"),
                  process_lvl = "proc_lvl")


# Create data frames for the different targets ----------------------------

source("analysis/scripts/indicators/createGraphDf.R")

species <- createGraphDf(elem_list, "species")
sp_type <- createGraphDf(elem_list, "sp_type")
site <- createGraphDf(elem_list, "site")
site_type <- createGraphDf(elem_list, "site_type")


# Create graphs for the targets -------------------------------------------

source("analysis/scripts/indicators/createGraph.R")

sp_gr <- createGraph(species)
sptpe_gr <- createGraph(sp_type)
sit_gr <- createGraph(site)
sitpe_gr <- createGraph(site_type)

render_graph(sp_gr)
render_graph(sptpe_gr)
render_graph(sit_gr)
render_graph(sitpe_gr)


# Combine graphs ----------------------------------------------------------

source("analysis/scripts/indicators/combineGraphs.R")

graph_comb <- combineGraphs(list(sp_gr, sit_gr, sptpe_gr, sitpe_gr))

render_graph(graph_comb)

export_graph(graph_comb, file_name = "analysis/scripts/indicators/all_graph_comb.png", file_type = "PNG")


# Find possible routes ----------------------------------------------------

ind_df <- site_type

# Find one route
source("analysis/scripts/indicators/createConditionGraph.R")

render_graph(createConditionGraph(ind_df, c("at", "across", NA), c("site_type", "sp_type", NA), output = "graph"))

# Find all routes
source("analysis/scripts/indicators/createAllConditionGraphs.R")

all_graphs <- createAllConditionGraphs(ind_df, output = "graph")

graph_comb <- combineGraphs(all_graphs)

render_graph(graph_comb)

export_graph(graph_comb, file_name = "analysis/scripts/indicators/site_graph_comb.png", file_type = "PNG")


# Export table ------------------------------------------------------------

inds_all <- bind_rows(
    species %>%
        mutate(taxon = NA, cond3 = cond2, cond2 = NA) %>%
        dplyr::select(names(sp_type)),
    sp_type %>%
        dplyr::select(names(sp_type)),
    site %>%
        mutate(geo = NA, cond3 = cond2, cond2 = NA) %>%
        dplyr::select(names(sp_type)),
    site_type %>%
        dplyr::select(names(sp_type))
)

write.csv(inds_all, file = "analysis/scripts/indicators/indtr_all.csv", row.names = FALSE)
