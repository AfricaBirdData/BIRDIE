#' Plot occupancy covariates
#'
#' @description Plots selected covariates against detection/non-detection data
#' @param occu_data An list with occupancy site and visit data. See \code{\link[BIRDIE]{ppl_create_site_visit}}
#' @param vars Character vector with the names of the variables to plot. They
#' must match those names in occu_data.
#' @param type The type of plot we want to display. Either "ind" which shows
#' individual plots for each variable, or "cross", which shows a scatterplot
#' with two variables. We need to pass on exactly two variables.
#'
#' @return A plot with detection/non-detection against selected variables
#' @export
#'
#' @examples
plotOccuVars <- function(occu_data, vars, type = "ind"){

    plotdata <- occu_data$visit %>%
        dplyr::left_join(occu_data$site %>%
                             dplyr::select(site, year, all_of(vars)),
                         by = c("site", "year")) %>%
        dplyr::select(obs, all_of(vars))

    if(type == "ind"){

        plotdata %>%
            tidyr::pivot_longer(cols = -obs, names_to = "vars", values_to = "value") %>%
            ggplot2::ggplot() +
            ggplot2::geom_jitter(aes(x = value, y = obs), height = 0.1, alpha = 0.3) +
            geom_smooth(aes(x = value, y = obs)) +
            ggplot2::ylab("Detection") +
            ggplot2::facet_wrap("vars", scales = "free")

    } else if(type == "cross"){

        plotdata %>%
            dplyr::mutate(obs = factor(obs)) %>%
            ggplot2::ggplot() +
            ggplot2::geom_point(aes_string(x = vars[1], y = vars[2], col = "obs"), alpha = 0.5)

    } else {
        error("Type must be either 'ind' or 'cross'")
    }

}
