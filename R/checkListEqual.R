#' Checks that all elements in a list are identical
#'
#' @param x A list
#'
#' @return TRUE if all elements in a list are identical, FALSE otherwise. If the function finds an element that is different it will print the elements that differ
#' @export
#'
#' @examples
#' checkListEqual(x = list(a = 1:5, b = 3:6))
#' checkListEqual(x = list(a = letters[7:9], b = c("g", "h", "i")))
checkListEqual <- function(x){

  for(i in 2:length(x)){

    out <- TRUE

    ref <- x[[i-1]]
    if(!identical(ref, x[[i]])){
      j <- i - 1
      print(paste(i, "differs from", j))
      out <- FALSE
      next
    }
  }

  return(out)

}
