print.lba.mle.fe <- function(x, 
                             digits = 3L, 
                             ...) {

  cat("\nCall:\n", 
      paste(deparse(x$call), 
            sep = "\n", 
            collapse = "\n"), 
      "\n")

  an <- round(x$A,
              digits)
  bn <- round(x$B, 
              digits)

  cat("\nIdentified mixing parameters:\n\n") 
  print.default(an,
                ...) 

  cat("\nIdentified latent budget:\n\n") 
  print.default(bn,
                ...) 

  cat("\nBudget proportions:\n") 
  pkk <- round(x$pk, 
               digits)
  rownames(pkk) <- c('')
  print.default(pkk, 
                ...) 
}
