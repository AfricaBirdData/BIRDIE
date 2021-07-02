#' Find index of next record matching criteria
#'
#' @description The function takes a vector and, for each element, finds the
#' index of the next record that matches certain criteria.
#' @param x A vector.
#' @param find_id A vector containing the elements that indices want to be
#' found of.
#'
#' @return A vector with the same length of x containg, for each element in x,
#' what is the index of the next record matching criteria passed on to find_id.
#' @export
#'
#' @examples
#' x <- rbinom(20, 1, 0.2)
#' findNextIndex(x, 1)
findNextIndex <- function(x, find_id){

    idx_w <- which(x %in% find_id)

    tgt_id <- vector("numeric", length = length(x))

    j <- 1

    for(i in seq_along(x)){

        if(j > length(idx_w)){
            tgt_id[i] <- NA
        } else if(i < idx_w[j]){
            tgt_id[i] <- idx_w[j]
        } else {
            j <- j + 1
            tgt_id[i] <- idx_w[j]
        }
    }

    return(tgt_id)
}
