goodnessfit.default <- function(object,...){

  stopifnot(inherits(object,'lba'))

  result <- goodnessfit(object, ...)

}
