#' Select occuR model
#'
#' @param forms A list with formulas to fit an occuR model see \link[occuR]{fit_occu}
#' @param type One of "visit" if the model selection is for the probability of
#' detection or "site" if model selection is for occupancy probability. An
#' overall selection to come...
#' @param visit_data Occupancy visit data produced by \link{prepDataOccuR} or with
#' the format prescribed by \link[occuR]{fit_occu}.
#' @param site_data Occupancy site data produced by \link{prepDataOccuR} or with
#' the format prescribed by \link[occuR]{fit_occu}.
#' @param mod_id Model identifier. It is recommended to include something that
#' allows us to identify which model and to what species.
#' species id.
#'
#' @return A dataframe with four columns: df - estimated degrees of freedom,
#' form - formula of the assessed model, AIC - AIC score of the model, mod -
#' model identifier.
#' @export
#'
#' @examples
selOccuRmod <- function(forms, type, visit_data, site_data, mod_id){

    if(type == "visit"){
        form <- forms[[1]]
    } else if(type == "site"){
        form <- forms[[2]]
    }

    fit <- fit_occu(forms = forms,
                    visit_data = visit_data,
                    site_data = site_data)

    out <- data.frame(df = dof.occuR(fit, each = FALSE),
                      form = Reduce(paste, deparse(form)),
                      AIC = AIC(fit),
                      mod = mod_id)

    return(out)

}
