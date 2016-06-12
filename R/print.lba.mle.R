print.lba.mle <- function(x, 
                          digits = 3L, 
                          ...) {

  cat("\nCall:\n", 
      paste(deparse(x$call), 
            sep = "\n", 
            collapse = "\n"), 
      "\n")

  aoi <- round(x$Aoi, 
               digits)
  boi <- round(x$Boi, 
               digits)

  cat("\nIdentified mixing parameters:\n\n") 
  print.default(aoi, 
                ...) 

  cat("\nIdentified latent budget:\n\n") 
  print.default(boi,
                ...)

  cat("\nBudget proportions:\n") 
  pkk <- round(x$pk, 
               digits)
  rownames(pkk) <- c('')
  print.default(pkk, 
                ...) 
}
