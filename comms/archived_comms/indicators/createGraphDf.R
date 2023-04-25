createGraphDf <- function(elem_list, target){

    groups <- c("taxon", "geo")
    group <- names(elem_list)[sapply(elem_list, function(x) target %in% x)]

    # if group geo, then remove global from filter options
    if(group == "geo"){
        elem_list$geo <- elem_list$geo[elem_list$geo != "global"]
    }

    elem_list$target <- target

    if(target %in% c("sp_type", "site_type")){

        elem_list$cond1 <- elem_list$cond
        elem_list$cond2 <- elem_list$cond
        elem_list$cond3 <- elem_list$cond

        elem_list <- elem_list[c("target", "cond1", groups[groups == group], "cond2", groups[groups != group], "cond3", "time", "indicator", "process_lvl")]

    } else if(target %in% c("species", "site")){

        elem_list$cond1 <- elem_list$cond
        elem_list$cond2 <- elem_list$cond

        elem_list <- elem_list[c("target", "cond1", groups[groups != group], "cond2", "time", "indicator", "process_lvl")]

    }

    graph_df <- expand.grid(elem_list, stringsAsFactors = FALSE)

    # eliminate non-necessary elements

    # More than one across condition
    n_across <- graph_df %>%
        select(starts_with("cond"))
    n_across <- rowSums(n_across == "across")

    graph_df <- graph_df %>%
        filter(!(n_across > 1))

    # Filters that match the target, with "across" condition
    if(sum(c("taxon", "geo") %in% names(graph_df)) > 1){

        graph_df <- graph_df %>%
            filter(!(target == taxon & cond1 == "across"),
                   !(target == geo & cond1 == "across"))

    }

    # Filters that don't match the target, with "at" condition (these should have their own diagrammes)
    if(sum(c("taxon", "geo") %in% names(graph_df)) > 1){

        if(group == "taxon"){
            graph_df <- graph_df %>%
                filter(!(target != taxon & cond1 == "at"))
        } else if(group == "geo"){
            graph_df <- graph_df %>%
                filter(!(target != geo & cond1 == "at"))
        }
    }

    # Remove across "global" conditions
    if(group == "taxon" && ncol(graph_df) == 7){
        graph_df <- graph_df %>%
            filter(!(cond1 == "across" & geo == "global"))
    } else if(group == "taxon" && ncol(graph_df) == 9){
        graph_df <- graph_df %>%
            filter(!(cond2 == "across" & geo == "global"))
    }

    # Set attributes to determine the possible conditions
    if(ncol(graph_df) == 7){
        attr(graph_df, "filters") <- names(graph_df)[c(3,5)]
    } else if(ncol(graph_df) == 9){
        attr(graph_df, "filters") <- names(graph_df)[c(3,5,7)]
    }

    return(graph_df)

}
