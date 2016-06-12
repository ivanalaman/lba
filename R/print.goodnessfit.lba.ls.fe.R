print.goodnessfit.lba.ls.fe <- function(x,
					digits=3L,
					...){

  cat("Residual sum of square:\n")

  res <- c(x$RSS1,
	   x$RSS)
  names(res) <- c("RSS baseline",
		  "RSS K budget")
  print.default(res,
		digits=digits,
		...)

  cat("\n")
}
