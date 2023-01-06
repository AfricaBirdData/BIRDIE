#' Create pipeline event log
#'
#' @description This function will create a customized log entry for a log file
#' that records the activity of the pipeline. More precisely, it will record the
#' time, species and outcome of data preparation, model fitting, diagnostics and
#' preparation.
#' @param config A list with pipeline configuration parameters.
#' See \code{\link{configPreambOccu}} and \code{\link{configPreambAbu}}.
#' @param logfile Optional. A character string with the path to the file that
#' should be updated. If NULL (default) a new log file will be created.
#' @param date_time Optional. A character string with the log date. If NULL (default)
#' the time given by \code{\link{Sys.time()}}, formatted as SAST will be used.
#' @param species SAFRING code of the species the log corresponds to.
#' @param model Type of model being processed. Either "occ" or "ssm".
#' @param data Data preparation status. Defaults to NA.
#' @param fit model fitting status. Defaults to NA.
#' @param diagnose Model diagnostics status. Defaults to NA.
#' @param summary Model summary status. Defaults to NA.
#' @param notes Any additional comments. Defaults to NA.
#' @param full_log Optionally, we can pass a full vector log, with all the above
#' elements. In this case, all other arguments but \code{config} and \code{logfile}
#' will be ignored.
#'
#' @return A .csv will be saved on the reports directory. The location of this
#' directory is configured in the configuration object passed as \code{config} at
#' \code{config$output_dir}.
#' @export
#'
#' @examples
createLog <- function(config, logfile = NULL, date_time = NULL, species = NULL,
                      model = NULL, data = NA, fit = NA, diagnose = NA,
                      summary = NA, notes = "", full_log = NULL){

    s <- Sys.time()
    attr(s,"tzone") <- "Africa/Johannesburg"
    today <- format(s, format = "%Y-%m-%d")

    if(is.null(date_time)){
        date_time <- format(s)
    }

    if(is.null(full_log)){
        new_log <- c(date_time, species, model, data, fit, diagnose, summary, notes)
    } else {
        new_log <- full_log
    }


    if(is.null(logfile)){
        logfile <- file.path(config$out_dir, "reports", paste0("pipe_log_", today,".csv"))
        log <- data.frame(date_time = character(),
                          species = numeric(),
                          model = character(),
                          data = numeric(),
                          fit = numeric(),
                          diagnose = numeric(),
                          summary = numeric(),
                          notes = character())
        utils::write.csv(log, logfile, row.names = FALSE)
        message(paste("log file created:", logfile))
    }

    log <- utils::read.csv(logfile)
    new_row <- nrow(log) + 1

    log[new_row, ] <- new_log

    utils::write.csv(log, logfile, row.names = FALSE)

}

