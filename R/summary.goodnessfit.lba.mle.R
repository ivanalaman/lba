summary.goodnessfit.lba.mle <- function(object,
                                        digits = 2L,
                                        ...){

  obj <- object

  stabasic <- round(matrix(c(obj[[4]],
                             obj[[2]],
                             obj[[8]],
                             obj[[3]],
                             obj[[1]],
                             obj[[7]]),
                           ncol=2),
                    digits)

  rownames(stabasic) <- c('Likelihood ratio test:',
                          'Degree of freedom:',
                          'Probability:')
  colnames(stabasic) <- c('LBM(K)','LBM(1)') 

  stachi <- round(matrix(c(obj[[6]],
                           obj[[10]], 
                           obj[[5]],
                           obj[[9]]),
                         ncol=2),
                  digits)

  rownames(stachi) <- c('Pearson Chi-square:',
                        'Probability:')
  colnames(stachi) <- c('LBM(K)','LBM(1)')  

  cms <- round(matrix(c(obj[[12]],
                        obj[[14]],
                        obj[[16]],
                        obj[[11]],
                        obj[[13]],
                        obj[[15]]),
                      ncol=2),
               digits)

  rownames(cms) <- c('AIC:',
                     'BIC:',
                     'CAIC:')
  colnames(cms) <- c('LBM(K)','LBM(1)')  

  ifi <- round(matrix(c(obj[[17]],
                        obj[[18]],
                        obj[[19]],
                        obj[[20]],
                        rep(NA,4)
                        ),
                      ncol=2),
               digits)

  rownames(ifi) <- c('Normed fit index:',
                     'Normed fit index modified:',
                     'Bollen index:',
                     'Tucker-Lewis index:')

  colnames(ifi) <- c('LBM(K)','LBM(1)')  

  otherm <- round(matrix(c(obj[[22]],
                           obj[[23]],
                           obj[[24]],
                           obj[[25]],
                           obj[[21]],
                           rep(NA,3)
                           ),
                         ncol=2),
                  digits)

  rownames(otherm) <- c('RSS:',
                        'Improvement:',
                        'Required per budget:',
                        'Required per defree of freedom:')
  colnames(otherm) <- c('LBM(K)','LBM(1)') 

  diss <- round(matrix(c(obj[[27]],
                         obj[[29]],
                         obj[[30]],
                         obj[[31]],
                         obj[[32]],
                         obj[[26]],
                         obj[[28]],
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

  madi <- round(matrix(c(obj[[34]],
                         obj[[35]],
                         obj[[36]],
                         obj[[37]],
                         obj[[33]],
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
              stachi=stachi,
              cms=cms,
              ifi=ifi,
              otherm=otherm,
              diss=diss,
              madi=madi)

  cat("LIKELIHOOD STATISTICS:\n\n")

  stab_out <- res$stabasic
  print.default(stab_out,
                na.print="-")

  cat("\nCHI-SQUARE STATISTICS:\n\n")

  stac_out <- res$stachi
  print.default(stac_out,
                na.print="-")

  cat("\nCRITERION FOR MODEL SELECTION:\n\n")

  cms_out  <- res$cms
  print.default(cms_out,
                na.print="-")

  cat("\nINCREMENTAL FIT INDICES:\n\n")

  ifi_out  <- res$ifi
  print.default(ifi_out,
                na.print="-")

  cat("\nBY RESIDUAL SUM OF SQUARE (RSS):\n\n")

  other_out<- res$otherm
  print.default(other_out,
                na.print="-")

  cat("\nBY DISSIMILARITY:\n\n")

  diss_out <- res$diss
  print.default(diss_out,
                na.print="-")

  cat("\nBY MEAN ANGULAR DEVIATION (MAD):\n\n")

  madi_out <- res$madi
  print.default(madi_out,
                na.print="-") 
  invisible(res)
}
