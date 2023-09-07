combineGraphs <- function(graphs, postprocess = FALSE){

    graph <- create_graph(attr_theme = "lr")

    n <- length(graphs)

    # Set initial nodes and edges
    maxnode <- 0
    new_node_df <- graphs[[1]]$nodes_df

    maxedge <- 0
    new_edge_df <- graphs[[1]]$edges_df

    # Define initial nodes
    ini_nodes <- c(target = 1, at = 2, across = 3)
    ini_nodes <- c(target = 1)

    # Loop through the rest of graphs and combine
    for(i in 2:n){

        maxnode <- max(new_node_df$id)

        new_node_df <- rbind(new_node_df,
                             graphs[[i]]$nodes_df %>%
                                 mutate(id = id + maxnode))



        maxedge <- max(new_edge_df$id)
        new_edge_df <- rbind(new_edge_df,
                             graphs[[i]]$edges_df %>%
                                 mutate(id = id + maxedge,
                                        to = to + maxnode,
                                        from = from + maxnode))

    }

    # Build new graph
    graph <- graph %>%
        add_node_df(new_node_df) %>%
        add_edge_df(new_edge_df)


    if(postprocess){

        # Post-processing ---------------------------------------------------------


        # Fix initial nodes

        new_node_df <- graph$nodes_df

        ini_ids <- list(target = new_node_df %>%
                            filter(type == "target", label == label[1]) %>%
                            pull(id),
                        at = new_node_df %>%
                            filter(type == "cond1", label == "at") %>%
                            pull(id),
                        across = new_node_df %>%
                            filter(type == "cond1", label == "across") %>%
                            pull(id))

        new_node_df <- new_node_df %>%
            mutate(id = if_else(type == "target" & label == label[1], 1L, id),
                   id = if_else(type == "cond1" & label == "at", 2L, id),
                   id = if_else(type == "cond1" & label == "across", 3L, id),
                   id = if_else(id == 3 & (type != "cond1" | label != "across"), max(id) + 1L, id))

        new_ini_ids <- new_node_df$id[1:3]

        # Fix initial edges
        new_edge_df <- graph$edges_df

        new_edge_df <- new_edge_df %>%
            mutate(to = if_else(to == 3, new_ini_ids[3], to),
                   from = if_else(from == 3, new_ini_ids[3], from))

        new_edge_df <- new_edge_df %>%
            mutate(to = case_when(to %in% ini_ids[[1]] ~ 1L,
                                  to %in% ini_ids[[2]] ~ 2L,
                                  to %in% ini_ids[[3]] ~ 3L,
                                  TRUE ~ to),
                   from = case_when(from %in% ini_ids[[1]] ~ 1L,
                                    from %in% ini_ids[[2]] ~ 2L,
                                    from %in% ini_ids[[3]] ~ 3L,
                                    TRUE ~ from))

        # Keep only one edge from 1 to (2,3)
        from1 <- new_edge_df %>%
            filter(from == 1) %>%
            pull(id)

        keep_from1 <- new_edge_df %>%
            filter(from == 1) %>%
            distinct(from, to, .keep_all = TRUE) %>%
            pull(id)

        new_edge_df <- new_edge_df %>%
            filter(!(id %in% setdiff(from1, keep_from1)))


        graph$nodes_df <- new_node_df
        graph$edges_df <- new_edge_df

    }

    return(graph)

}
