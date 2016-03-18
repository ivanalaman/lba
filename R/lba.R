lba <- function (obj, 
                 ...) UseMethod('lba')

print.lba <- function(x,
                      digits = max(3L, getOption("digits") - 3L),
                      ...)
{
 cat("\nCall:\n",
     paste(deparse(x$call),
           sep = "\n",
           collapse = "\n"),
     "\n\n",
     sep = "")

 cat("Mixing parameters:\n",
 ifelse(class(x)[1] == 'lba.ls' | class(x)[1] == 'lba.mle',
        print.default(format(x$Aoi,
                             digits = digits),
                      print.gap = 2L,
                      quote = FALSE),
        print.default(format(x$A,
                             digits = digits),
                      print.gap = 2L,
                      quote = FALSE)),
     "\n\n",
     sep = "")

 cat("Latent budgets:\n")
 ifelse(class(x)[1] == 'lba.ls' | class(x)[1] == 'lba.mle',
        print.default(format(x$Boi,
                             digits = digits),
                      print.gap = 2L,
                      quote = FALSE),
        print.default(format(x$B,
                             digits = digits),
                      print.gap = 2L,
                      quote = FALSE))
  
 cat("\n")
 invisible(x)

}
