#' Predict from occuR model fit
#'
#' @inheritParams ppl_run_pipe_distr
#'
#' @return
#' @export
#'
#' @examples
ppl_predict_occur <- function(sp_code, year, config, ...){

    varargs <- list(...)

    occuRdata <- ppl_prep_occur_data(sp_code, year, config, ...)


    # Prepare prediction data -------------------------------------------------

    # Load site data for prediction
    if(year < 2020){
        datafile <- file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")
    } else {
        yy <- substring(as.character(config$year), 3, 4)
        ff <- paste0("site_dat_sa_gee_", yy, ".rds")
        datafile <- file.path(config$data_dir, ff)
    }

    # Load data, save geometry and subset years
    predsites <- readRDS(datafile)

    gm <- predsites %>%
        dplyr::select(Pentad = Name, geometry)

    predsites <- predsites %>%
        sf::st_drop_geometry() %>%
        dplyr::select(Pentad = Name, lon, lat, watocc_ever, dist_coast, ends_with(match = as.character(config$years))) %>%
        tidyr::drop_na()   # I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE


    # Predict from model ------------------------------------------------------

    # Load fit
    fit <- with(config,
                readRDS(file.path(config$fit_dir, sp_code, paste0("occur_fit_", years_ch, "_", sp_code, ".rds"))))

    # Prepare prediction data
    pred_data <- prepPredictDataOccuR(occuRdata = occuRdata,
                                      pred_sites = predsites,
                                      years = config$years,
                                      scaling = TRUE)

    print(paste0("Predicting from model at ", Sys.time()))

    tryCatch({

        # Predict
        pred_occu <- predict(fit, occuRdata$visit, pred_data, nboot = 1000)

    }, error = function(e){
        sink(file.path(config$fit_dir, sp_code, paste0("failed_pred_", sp_code,".txt")))
        print(e)
        sink()}) # TryCatch predict

    return(list(pred = pred_occu,
                occuRdata = occuRdata,
                pred_data = pred_data,
                pred_sites = predsites %>%
                    dplyr::left_join(gm, by = "Pentad") %>%
                    sf::st_as_sf()))

}
