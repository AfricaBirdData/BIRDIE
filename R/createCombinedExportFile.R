#' Create a combined file for export to data mart
#'
#' @param config A list with pipeline configuration parameters.
#' See \link{configPreambJAGS} and \link{configPreambOccuR}
#' @param type A character string with three options: "abu", for abundance estimates files,
#' "dst" for occupancy estimates files, or "indtr_dst" for distribution indicator files.
#'
#' @return It will create a file in config$out_dir combining all files of the selected
#' type of all species and years in config$years.
#' @export
#'
#' @examples
createCombinedExportFile <- function(config, type = c("abu", "dst", "indtr_dst")){

    if(type == "abu"){

        abu_out <- data.frame()

        for(i in 1:length(config$species)){

            sp_code <- config$species[i]

            abu_file <- file.path(config$out_dir, sp_code, paste0("ssm_pred_", config$years_ch, "_", sp_code, "_all.csv"))

            abu_out <- rbind(abu_out, read.csv(abu_file))
            if(file.exists(abu_file)){
                abu_out <- rbind(abu_out, read.csv(abu_file))
            } else {
                abu_out <- abu_out
            }

        }

        utils::write.csv(abu_out, file.path(config$out_dir, "export", paste0("ssm_pred_", config$years_ch, "_all_all.csv")),
                         row.names = FALSE)
    }

    if(type == "dst"){

        dst_out <- data.frame()

        for(i in 1:length(config$species)){

            sp_code <- config$species[i]

            for(y in 1:length(config$years)){

                yr <- substring(as.character(config$years[y]), 3, 4)
                dst_file <- file.path(config$out_dir, sp_code, paste0("occur_pred_", yr, "_", sp_code, ".csv"))

                if(file.exists(dst_file)){
                    dst_out <- rbind(dst_out, read.csv(dst_file))
                } else {
                    dst_out <- dst_out
                }

            }

        }

        utils::write.csv(dst_out, file.path(config$out_dir, "export", paste0("occur_pred_", config$years_ch, "_all.csv")),
                         row.names = FALSE)
    }

    if(type == "indtr_dst"){

        dst_out <- data.frame()

        for(i in 1:length(config$species)){

            sp_code <- config$species[i]

            dst_file <- file.path(config$out_dir, sp_code, paste0("indtr_dst_", sp_code, "_", config$year, ".csv"))

            if(file.exists(dst_file)){
                dst_out <- rbind(dst_out, read.csv(dst_file))
            } else {
                dst_out <- dst_out
            }

        }

        utils::write.csv(dst_out, file.path(config$out_dir, "export", paste0("indtr_dst_all_", config$years_ch, ".csv")),
                         row.names = FALSE)
    }


}

