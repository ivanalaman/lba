plotlba.lba.2d <- function(x,
                           axis.labels    = NULL,
                           labels.points  = NULL,
                           col.points     = par('fg'),
                           pch.budget     = par('pch'),
                           col.budget     = par('fg'),
                           lty.budget     = par('lty'),
                           lwd.budget     = par('lwd'),
                           colline.budget = par('fg'), 
                           args.legend    = NULL, 
                           with.ml        = c("mix","lat"),
                           ...)
{

  switch(match.arg(with.ml),
         mix = {
           if(inherits(x, 'lba.ls.fe') | inherits(x, 'lba.mle.fe') | inherits(x, 'lba.ls.logit') | inherits(x, 'lba.mle.logit')) {
             nrowss <- dim(x[[4]])[1]
             alfas <- x[[4]]
             alfas <- alfas/rowSums(alfas)
             pk <- x[[7]]
           } else {
             nrowss <- dim(x[[6]])[1]
             alfas <- x[[6]]
             alfas <- alfas/rowSums(alfas)
             pk <- x[[9]] 
           }
         },
         lat = {
           if(inherits(x, 'lba.ls.fe') | inherits(x, 'lba.mle.fe') | inherits(x, 'lba.ls.logit') | inherits(x, 'lba.mle.logit')) { 
             nrowss <- dim(x[[5]])[1]

             alfas <- x[[6]] 
             alfas <- alfas/rowSums(alfas) 

             pk <- x[[7]]

           } else {
             nrowss <- dim(x[[7]])[1]

             alfas <- x[[8]] 
             alfas <- alfas/rowSums(alfas) 

             pk <- x[[9]] 

           }
         }
         )

  rnames <- rownames(alfas) 
  K <- ncol(alfas)  

  if(K==1 | K>3) stop('Only K between 2 and 3 budgets') 

  if(is.null(labels.points)) labels.points <- as.character(1:nrowss) 

  invalfas <- alfas[,rev(order(colnames(alfas)))]
  invpk <- matrix(pk[,rev(order(colnames(pk)))],
                  ncol = 3)

  triax.plot(x            = invalfas,
             show.legend  = FALSE,
             label.points = TRUE,
             point.labels = '',
             pch          = '',
             axis.labels  = axis.labels,
             ...)

  triax.points(x = invpk,
               show.legend = FALSE,
               col.symbols = col.budget,
               pch         = pch.budget)

  xpos <- 1 - (invalfas[,1] + invalfas[,3] * 0.5)
  ypos <- invalfas[,3] * sin(pi/3)

  thigmophobe.labels(x = xpos,
                     y = ypos,
                     labels = labels.points,
                     col = col.points,
                     ...)

  segments(x0  = 0, 
           y0  = 0, 
           x1  = 1-(pk[,3] + pk[,1] * 0.5), 
           y1  = pk[,1] * sin(pi/3),
           lty = lty.budget,
           lwd = lwd.budget,
           col = colline.budget)

  segments(x0  = 1,
           y0  = 0, 
           x1  = 1-(pk[,3] + pk[,1] * 0.5), 
           y1  = pk[,1] * sin(pi/3),
           lty = lty.budget,
           lwd = lwd.budget,
           col = colline.budget)

  segments(x0  = 1-(pk[,3] + pk[,1] * 0.5), 
           y0  = pk[,1] * sin(pi/3), 
           x1  = 0.5, 
           y1  = 0.865,
           lty = lty.budget,
           lwd = lwd.budget,
           col =colline.budget)

  if(is.null(args.legend)){

    thelabels <- paste(as.character(1:nrowss),
                       ' ',
                       rnames,
                       sep="")

    args.2Kl <- list(x        = 'topleft',
                     legend   = thelabels,
                     text.col = col.points)

    do.call('legend',
            args.2Kl)

  } else {

    thelabels <- paste(as.character(1:nrowss),
                       ' ',
                       rnames,
                       sep="")

    args.2Kl <- list(x        = 'topleft',
                     legend   = thelabels,
                     text.col = col.points)

    args.2Kl[names(args.legend)] <- args.legend     

    do.call('legend',
            args.2Kl) 

  }    
}
