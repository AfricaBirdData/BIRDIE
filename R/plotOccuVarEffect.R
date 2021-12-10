#' Plot effect of occupancy covariates
#'
#' @description Plots the effect of selected covariates
#' @param fit An occuR model fit.
#' @param new_occu_data An list with occupancy site and visit data to predict
#' from. See \code{\link[BIRDIE]{prepDataOccuR}}
#' @param var Character string with the name of the variables to plot. It
#' must match those names in new_occu_data.
#' @param type A character string of either "occu" if the variable is plotted
#' against occupancy probabilities or "detc" if the variable is plotted against
#' detection probabilities
#' @param nboot Number of bootstrap samples to compute confidence intervals.
#'
#' @return A plot of the effect of selected variables on occupancy or detection.
#' probabilities.
#' @export
#'
#' @examples
plotOccuVarEffect <- function(fit, new_occu_data, var, type = "occu", nboot = 1000){

    pred_temp <- predict(obj = fit,
                         visit_data = new_occu_data$visit,
                         site_data = new_occu_data$site,
                         nboot = nboot)

    if(type == "occu" && nboot > 0){

        ci <- apply(pred_temp$psiboot, 2, quantile, prob = c(0.025, 0.975))

        data.frame(est = pred_temp$psi[,1],
                   lb = ci[1,],
                   ub = ci[2,],
                   var = dplyr::pull(new_occu_data$site, var)) %>%
            tidyr::pivot_longer(cols = -var,
                                names_to = "bound",
                                values_to = "value") %>%
            ggplot() +
            geom_line(aes(x = var, y = value, linetype = bound)) +
            ylab("Occupancy probability") + xlab(var)

    } else if(type == "occu" && nboot == 0){

        data.frame(est = pred_temp$psi,
                   var = dplyr::pull(new_occu_data$site, var)) %>%
            ggplot() +
            geom_line(aes(x = var, y = est)) +
            ylab("Occupancy probability") + xlab(var)

    } else if(type == "detc"){

        plotdata %>%
            dplyr::mutate(obs = factor(obs)) %>%
            ggplot2::ggplot() +
            ggplot2::geom_point(aes_string(x = vars[1], y = vars[2], col = "obs"), alpha = 0.5)

    } else {
        error("Type must be either 'occu' or 'detc'")
    }

}
