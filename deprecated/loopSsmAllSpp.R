#' Loop a state-space model through all species
#'
#' @description A wrapper around fitCwacSsm that runs a model for each species
#' and produces specific output.
#' @param counts A data frame with at least four columns: "counts" - an integer
#' column corresponding to the counts of a year and season, "season_id" - an
#' integer column that identifies the season (1 for summer and 2 for winter),
#' "spp" the species codes and "Loc_Code" the code for the site the counts where
#' taken at.
#' @param mod_file A character string corresponding to the directory where the
#' JAGS model lives at. If not provided, writeJagsModelFile is called to create
#' a "default" model.
#' @param data_outdir Directory where data outputs should be saved
#' @param plot_outdir Directory where plot outputs should be saved
#' @param ... Other arguments passed on to fitCwacSsm.
#'
#' @return The function returns the results from fitting a state-space model to
#' each of the species in two formats: i) a dataframe with summer and winter
#' states, as well as, trends and winter proportions, and ii) plots showing the
#' evolution of these variables.
#' @export
#'
#' @examples
#' #counts <- barberspan
#' #counts <- counts[counts$spp %in% c(6, 41),]
#' #loopSsmAllSpp(counts,
#' #              data_outdir = "mydatadir/",
#' #              plot_outdir = "myoutdir/",
#' #              param = c("beta", "lambda", "sig.zeta", "sig.w",
#' #               "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))
loopSsmAllSpp <- function(counts, mod_file = NULL, data_outdir, plot_outdir, ...){

    if(is.null(mod_file)){
        modpath <- BIRDIE::writeJagsModelFile()
        warning("mod_file not provided, running writeJagsModelFile, see ?writeJagsModelFile")
    } else {
        modpath <- mod_file
    }

    # Identify site
    site_id <- unique(counts$Loc_Code)

    if(length(site_id) != 1){
        stop("Either location code is missing or there are multiple sites.")
    }

    # Get species list
    spp <- unique(counts$spp)

    for(i in seq_along(spp)){

        sp <- spp[i]

        # Prepare data to fit an SSM
        ssmcounts <- BIRDIE::prepSsmData(counts, species = sp)

        # Fit 2-season dynamic trend model
        fit_dyn <- BIRDIE::fitCwacSsm(ssmcounts, mod_file = modpath, ...)

        # Plot
        p <- BIRDIE::plotSsm2ss(fit = fit_dyn, ssm_counts = ssmcounts, dyn = TRUE)

        # Extract plot data
        plotdata <- p$data

        # Prepare state data
        stt_df <- plotdata[[1]]

        stt_dfw <- cbind(stt_df %>%
                             dplyr::filter(season == 1) %>%
                             dplyr::select(-season),
                         stt_df %>%
                             dplyr::filter(season == 2) %>%
                             dplyr::select(-c(year, season)))

        names(stt_dfw) <- c("summer.logest", "summer.logest.ci.lower", "summer.logest.ci.upper",
                            "year", "log.summer.count",
                            "winter.logest", "winter.logest.ci.lower", "winter.logest.ci.upper",
                            "log.winter.count")

        # Prepare trend data
        trd_df <- plotdata[[2]]

        names(trd_df) <- c("slope.est", "slope.ci.lower", "slope.ci.upper",
                           "winter.prop.est", "winter.prop.ci.lower", "winter.prop.ci.upper",
                           "year")

        # Combine
        out_df <- cbind(stt_dfw, dplyr::select(trd_df, -year))

        # Add missing columns
        out_df <- out_df %>%
            dplyr::mutate(species = sp,
                          summer.count = round(exp(log.summer.count)),
                          winter.count = round(exp(log.winter.count)))

        # Order columns
        out_df <- out_df %>%
            dplyr::select(species, year, summer.count, winter.count,
                          log.summer.count, log.winter.count,
                          summer.logest, winter.logest, winter.logest.ci.lower,
                          winter.logest.ci.upper,
                          summer.logest.ci.lower, summer.logest.ci.upper,
                          slope.est, slope.ci.lower, slope.ci.upper,
                          winter.prop.est, winter.prop.ci.lower,
                          winter.prop.ci.upper)

        # Export sample
        datadir <- paste0(data_outdir, sp, "/ssm_dat_", site_id, "_", sp, ".csv")
        utils::write.csv(out_df, datadir, row.names = FALSE)

        plotdir <- paste0(plot_outdir, sp, "/ssm_plot_", site_id, "_", sp, ".png")
        ggplot2::ggsave(plotdir, plot = p$plot, width = 6.4, height = 5.2)
    }

    # Remove temporary model file
    if(is.null(mod_file)){
        unlink(modpath)
    }
}
