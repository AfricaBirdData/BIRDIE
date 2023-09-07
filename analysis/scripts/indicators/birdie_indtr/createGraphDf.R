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

        elem_list <- elem_list[c("target", "indicator", "cond1", groups[groups == group], "cond2", groups[groups != group], "cond3", "time")]

    } else if(target %in% c("species", "site")){

        elem_list$cond1 <- elem_list$cond
        elem_list$cond2 <- elem_list$cond

        elem_list <- elem_list[c("target", "indicator", "cond1", groups[groups != group], "cond2", "time")]

    }

    graph_df <- expand.grid(elem_list, stringsAsFactors = FALSE)


    # Eliminate non-necessary elements ----------------------------------------

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

    # Filters that don't match the target, with "at" condition (these should have their own diagrams)
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
    if(group == "taxon" && ncol(graph_df) == 6){
        graph_df <- graph_df %>%
            filter(!(cond1 == "across" & geo == "global"))
    } else if(group == "taxon" && ncol(graph_df) == 8){
        graph_df <- graph_df %>%
            filter(!(cond2 == "across" & geo == "global"))
    }

    # Remove distribution indicators for site and site_type
    if(target %in% c("site", "site_type")){
        graph_df <- graph_df %>%
            filter(indicator != "distr")
    }

    # Remove ancillary indicators for species and sp_type
    if(target %in% c("species", "sp_type")){
        graph_df <- graph_df %>%
            filter(indicator != "ancll")
    }

    # Remove diversity indicators for species
    if(target == "species"){
        graph_df <- graph_df %>%
            filter(indicator != "div")
    }

    if(target != "species"){
        graph_df <- graph_df %>%
            filter(!(indicator == "div" & taxon == "species"))
    }

    # Set attributes to determine the possible conditions
    if(ncol(graph_df) == 6){
        attr(graph_df, "filters") <- names(graph_df)[c(4,6)]
    } else if(ncol(graph_df) == 8){
        attr(graph_df, "filters") <- names(graph_df)[c(4,6,8)]
    }

    return(graph_df)

}
