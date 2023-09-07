#' Predict from occuR model fit
#'
#' @inheritParams ppl_run_pipe_distr
#'
#' @return
#' @export
#'
#' @examples
predictOccuR <- function(sp_code, year, config, ...){

    varargs <- list(...)

    occuRdata <- prepOccuRData(sp_code, year, config, ...)

    # Stop if there are too few detections
    if(is.numeric(occuRdata) && occuRdata %in% c(1, 2)){
        return(1)
    }


    # Prepare prediction data -------------------------------------------------

    # Load site data for prediction
    if(year < 2020){
        sitefile <- file.path(config$out_dir, "site_dat_sa_gee_08_19.csv")
    } else {
        sitefile <- file.path(config$out_dir, paste0("site_dat_sa_gee_", config$years_ch, ".csv"))
    }

    # Load data, save geometry and subset years
    predsites <- read.csv(sitefile) %>%
        dplyr::select(Pentad = Name, lon, lat, watocc_ever, dist_coast, dplyr::ends_with(match = as.character(config$years))) %>%
        tidyr::drop_na()   # I'M REMOVING SITES WITH NA DATA! MAKE SURE THIS MAKES SENSE


    # Predict from model ------------------------------------------------------

    # Load fit
    fit <- with(config,
                readRDS(file.path(out_dir, sp_code, paste0("occur_fit_", years_ch, "_", sp_code, ".rds"))))

    # Prepare prediction data
    pred_data <- prepPredictOccuRData(occuRdata = occuRdata,
                                      pred_sites = predsites,
                                      years = config$years,
                                      scaling = TRUE)

    print(paste0("Predicting from model at ", Sys.time()))

    tryCatch({

        # Predict
        pred_occu <- predict(fit, occuRdata$visit, pred_data, nboot = 1000)

    }, error = function(e){
        sink(file.path(config$out_dir, sp_code, paste0("failed_pred_", sp_code,"_", year, ".txt")))
        print(e)
        sink()}) # TryCatch predict

    return(list(pred = pred_occu,
                occuRdata = occuRdata,
                pred_data = pred_data,
                pred_sites = predsites)
           )

}
