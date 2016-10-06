plotlba.lba.1d <- function(x, 
                           height.line    = NULL,
                           xlab           = NULL,
                           ylab           = NULL,
                           ylim           = NULL,
                           args.legend    = NULL,
                           labels.points  = NULL,
                           col.points     = par('col'),
                           col.lines      = par('col'),
                           lty.lines      = par('lty'),
                           lwd.lines      = par('lwd'),
                           pch.budget     = par('pch'),
                           col.budget     = par('fg'),
                           lty.budget     = par('lty'),
                           lwd.budget     = par('lwd'),
                           colline.budget = NULL, 
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

  if(K==1) stop('Only K between 2 and 3 budgets') 

  if(is.null(height.line)) height.line <- seq(0,0.8,length=nrowss)

  if(is.null(xlab)) xlab <- ''

  if(is.null(ylab)) ylab <- ''

  if(is.null(ylim)) ylim <- c(0,1)

  if(is.null(labels.points)) labels.points <- as.character(1:nrowss)

  if(is.null(colline.budget)) colline.budget <- rgb(047,079,079,maxColorValue=255)

  plot(x    = alfas[,1],
       y    = height.line,
       xlab = xlab,
       ylab = ylab,
       ylim = ylim,
       type = 'h',
       axes = F,
       col  = col.lines,
       lty  = lty.lines,
       lwd  = lwd.lines,
       ...) 
  
  axis(1,
       ...)
  
  text(x      = alfas[,1],
       y      = height.line+0.02,
       labels = labels.points,
       col    = col.points,
       ...)

  points(pk[1,1],
         1,
         pch = pch.budget,
         col = col.budget,
         ...)

  segments(x0 = pk[1,1],
           y0 = 0,
           x1 = pk[1,1],
           y1 = 1,
           lty = lty.budget,
           lwd = lwd.budget,
           col = colline.budget)

  if(is.null(args.legend)){

    thelabels <- paste(as.character(1:nrowss),
                       ' ',
                       rnames,
                       sep="")

    args.2Kl <- list(x        = -0.2,
                     y        = 1.2,
                     legend   = thelabels,
                     text.col = col.points,
                     xpd = T)

    do.call('legend',
            args.2Kl)

  } else {

    thelabels <- paste(as.character(1:nrowss),
                       ' ',
                       rnames,
                       sep="")

    args.2Kl <- list(x        = -0.2,
                     y        = 1.2,
                     legend   = thelabels,
                     text.col = col.points,
                     xpd=T)

    args.2Kl[names(args.legend)] <- args.legend     

    do.call('legend',
            args.2Kl) 

  }       
} 
