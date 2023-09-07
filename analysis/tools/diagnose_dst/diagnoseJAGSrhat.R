diagnoseJAGSrhat <- function(fit, param = "all"){

    if(param == "all"){
        rhats <- lapply(fit$Rhat, function(x) sum(abs(x - 1) > 0.1)) %>%
            as.data.frame()
    }

    if(param != "all"){

        if(length(param) != 1){
            stop("Currently, we can only assess one parameter at a time")
        }

        # Extract rhats
        rhats <- eval(parse(text = paste0("fit$Rhat$", param)))

        if(length(dim(rhats)) == 1){

            # Plot
            rhatPlot <- data.frame(rhat = rhats) %>%
                dplyr::mutate(year = dplyr::row_number(),
                              convergence = factor(ifelse(abs(rhats - 1) > 0.1, 0, 1))) %>%
                ggplot2::ggplot() +
                ggplot2::geom_point(ggplot2::aes(x = year, y = rhat,
                                        col = convergence))

            # Find specific sites with non-compliant Rhats
            year_ids_rhat <- data.frame(rhat = rhats) %>%
                dplyr::mutate(year = dplyr::row_number(),
                              convergence = factor(ifelse(abs(rhats - 1) > 0.1, 0, 1))) %>%
                dplyr::filter(convergence == 0) %>%
                dplyr::pull(year) %>%
                unique()

            message(paste("Non-convergence on year", paste(year_ids_rhat, collapse = " ")))

        } else {
            colnames(rhats) <- paste0("year_", seq_len(ncol(rhats)))
            rhats <- rhats %>%
                as.data.frame() %>%
                dplyr::mutate(site = dplyr::row_number())

            # Plot
            rhatPlot <- rhats %>%
                tidyr::pivot_longer(cols = -site, names_to = "year", values_to = "rhat") %>%
                dplyr::mutate(year = gsub("year_", "", year)) %>%
                dplyr::mutate(year = factor(year, levels = unique(year)),
                              site = factor(site, levels = unique(site))) %>%
                ggplot2::ggplot() +
                ggplot2::geom_raster(ggplot2::aes(x = site, y = year, fill = abs(rhat-1)))

            # Find specific sites with non-compliant Rhats
            site_ids_rhat <- rhats %>%
                tidyr::pivot_longer(cols = -site, names_to = "year", values_to = "rhat") %>%
                dplyr::mutate(year = gsub("year_", "", year)) %>%
                dplyr::mutate(year = factor(year, levels = unique(year)),
                              site = factor(site, levels = unique(site))) %>%
                dplyr::filter(abs(rhat - 1) > 0.1) %>%
                dplyr::pull(site) %>%
                unique() %>%
                as.numeric()

            message(paste("Non-convergence on sites", paste(site_ids_rhat, collapse = " ")))
        }
    }

    print(plot(rhatPlot))

}
