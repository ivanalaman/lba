summary.goodnessfit.lba.ls.logit <- function(object,
                                             digits = 2L,
                                             ...){

  obj <- object

  stabasic <- round(matrix(c(obj[[2]],
                             obj[[1]]),
                           ncol=2),
                    digits)

  rownames(stabasic) <- 'Degree of freedom:'
  colnames(stabasic) <- c('LBM(K)','LBM(1)') 

  otherm <- round(matrix(c(obj[[4]],
                           obj[[5]],
                           obj[[6]],
                           obj[[7]],
                           obj[[3]],
                           rep(NA,3)
                           ),
                         ncol=2),
                  digits)

  rownames(otherm) <- c('RSS:',
                        'Improvement:',
                        'Required per budget:',
                        'Required per defree of freedom:')
  colnames(otherm) <- c('LBM(K)','LBM(1)')  

  diss <- round(matrix(c(obj[[9]],
                         obj[[11]],
                         obj[[12]],
                         obj[[13]],
                         obj[[14]],
                         obj[[8]],
                         obj[[10]],
                         rep(NA,3)
                         ),
                       ncol=2),
                digits)

  rownames(diss) <- c('Index of dissimilarity:',
                      'Prop. correctly classf. data:',
                      'Improvement:',
                      'Required per budget:',
                      'Required per defree of freedom:')    

  colnames(diss) <- c('LBM(K)','LBM(1)')  

  madi <- round(matrix(c(obj[[16]],
                         obj[[17]],
                         obj[[18]],
                         obj[[19]],
                         obj[[15]],
                         rep(NA,3)
                         ),
                       ncol=2),
                digits)

  rownames(madi) <- c('MAD:',
                      'Improvement:',
                      'Required per budget:',
                      'Required per defree of freedom:')    
  colnames(madi) <- c('LBM(K)','LBM(1)')  

  res <- list(stabasic=stabasic,
              otherm=otherm,
              diss=diss,
              madi=madi)

  cat("BASIC STATISTICS:\n\n")

  stab_out <- res$stabasic
  print.default(stab_out)

  cat("\nBY RESIDUAL SUM OF SQUARE (RSS):\n\n") 

  outherm_out <- res$otherm
  print.default(outherm_out,
                na.print="-")

  cat("\nBY DISSIMILARITY:\n\n")

  diss_out <- res$diss
  print.default(diss_out,
                na.print="-")

  cat("\nBY MEAN ANGULAR DEVIATION (MAD):\n\n")

  mat_out <- res$madi   
  print.default(mat_out,
                na.print="-") 

}
