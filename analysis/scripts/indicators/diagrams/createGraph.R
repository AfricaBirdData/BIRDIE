createGraph <- function(df){

    # Set up general graph options --------------------------------------------

    graph <- create_graph(attr_theme = "lr")


    # Add node data frames -------------------------------------------------

    # Node colours
    col_sel <- viridis::turbo(length(names(df)))

    # Font colours
    n_names <- length(names(df))
    font_col_sel <- c("white", rep("black", n_names - 2), "white")

    for(i in seq_along(names(df))){

        cat <- names(df)[i]

        new_nodes <- create_node_df(
            n = n_distinct(df[,cat]),
            type = cat,
            label = unique(df[,cat]),
            fontcolor = font_col_sel[i],
            fillcolor = col_sel[i],
            color = "gray70",
            penwidth = 1,
            width = 0.7
        )

        graph <- graph %>%
            add_node_df(new_nodes)

    }


    # Add edges data frames ------------------------------------------------

    get_node_df(graph)

    for(i in seq_along(names(df)[-1])){

        cat1 <- names(df)[i]
        cat2 <- names(df)[i+1]
        lvs <- unique(df[,cat1])
        lvs <- lvs[!is.na(lvs)]

        for(j in seq_along(lvs)){

            d <- lvs[j]

            from <- df %>%
                filter(df[, cat1] == d) %>%
                distinct(.[, cat1], .[, cat2]) %>%
                rename(label = 1) %>%
                left_join(get_node_df(graph) %>%
                              filter(type == cat1),
                          by = "label") %>%
                pull(id)

            to <- df %>%
                filter(df[, cat1] == d) %>%
                distinct(.[, cat1], .[, cat2]) %>%
                rename(label = 2) %>%
                left_join(get_node_df(graph) %>%
                              filter(type == cat2),
                          by = "label") %>%
                pull(id)

            new_edges <- create_edge_df(
                from = from,
                to = to,
                rel = "toward",
                color = "gray80",
                penwidth = 1)

            graph <- graph %>%
                add_edge_df(new_edges)

        }
    }

    return(graph)

}
