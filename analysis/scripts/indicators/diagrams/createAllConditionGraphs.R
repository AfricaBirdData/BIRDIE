createAllConditionGraphs <- function(df, output = "graph"){

    filts <- attr(df, "filters")
    nfilts <- length(filts)

    # The last filter is redundant, so we can remove
    filts <- filts[-nfilts]

    # Find the indexes in the df and the predecing conditions
    idx <- which(names(df) %in% filts)
    idx <- sort(c(idx, idx - 1))

    # Find possible combinations
    n_comb <- n_distinct(df[,idx])
    comb <- distinct(df[,idx])

    graph_list <- vector("list", length = n_comb)

    for(i in 1:n_comb){
        graph_list[[i]] <- createConditionGraph(df,
                                                comb %>%
                                                    select(starts_with("cond")) %>%
                                                    slice(i) %>%
                                                    as.vector() %>%
                                                    append(NA),
                                                comb %>%
                                                    select(!starts_with("cond")) %>%
                                                    slice(i) %>%
                                                    as.vector() %>%
                                                    append(NA),
                                                output = output)
    }

    return(graph_list)

}
