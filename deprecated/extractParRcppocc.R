#' Extract samples from a parameter in Rcppocc multi-chain fit
#'
#' @param fit Multichain Rcppocc fit obtained via prllOccSPATlogit
#' @param par Parameter we want to extract samples from
#'
#' @return A dataframe with the samples of the parameter of interest with a
#' column indicating the chain the samples correspond to.
#' @export
#'
#' @examples
extractParRcppocc <- function(fit, pars){

    smp <- fit[[pars]]

    varnames <- paste0(pars, seq_len(nrow(smp)))
    nchain <- dim(smp)[3]

    smp <- array(smp, dim = c(dim(smp)[1], dim(smp)[2]*nchain))

    smp <- as.data.frame(t(smp))
    names(smp) <- varnames

    smp$chain <- rep(seq_len(nchain), each = nrow(smp)/nchain)

    return(smp)

    }

