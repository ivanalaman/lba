print.goodnessfit.lba.mle.logit <- function(x,
					    digits=3L,
					    ...){

  cat("Likelihood ratio statistic:\n")

  mat <- matrix(c(x$G2,
		  x$proG,
		  x$G2b,
		  x$proG1),
		ncol=2)
  rownames(mat) <- c('G2 value',
		     'P-value') 
  colnames(mat) <- c('K budget',
		     'Baseline')

  print.default(mat,
		digits=digits,
		...) 
}
