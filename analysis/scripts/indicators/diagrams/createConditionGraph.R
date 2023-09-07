createConditionGraph <- function(df, cond, filter, output = "graph"){

    filts <- attr(df, "filters")
    nfilts <- length(filts)

    if(nfilts != length(cond) || nfilts != length(filter)){
        stop(paste("You need to provide", nfilts, "conditions and filters. Set to NA those you don't want to use."))
    }

    print(paste(sum(!is.na(cond)), "of", nfilts, "conditions available set"))

    for(i in which(!is.na(cond))){

        df <- df[df[,paste0("cond", i)] == cond[i],]
        df <- df[df[,filts[i]] == filter[i],]

    }

    if(nrow(df) == 0){
        stop("There is nothing with those specifications")
    }

    # if(cond == "across"){
    #     new_df <- new_df %>%
    #         filter(cond2 != cond)
    # }

    if(output == "df"){

        out <- df

    } else if(output == "graph"){

        out <- createGraph(df)

    }


    return(out)

}
