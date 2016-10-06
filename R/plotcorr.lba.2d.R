plotcorr.lba.2d <- function(x,
                         dim            = c(1,2), #only K = 3
                         xlim           = NULL,
                         ylim           = NULL,
                         xlab           = NULL,
                         ylab           = NULL,
                         args.legend    = NULL,
                         col.points     = NULL,
                         labels.points  = NULL,
                         pch.points     = NULL,
                         pos.points     = NULL,
                         labels.budget  = NULL,
                         pch.budget     = NULL,
                         pos.budget     = NULL,
                         cex.budget     = NULL,
                         col.budget     = NULL,
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
             alfas <- x[[5]] 
             pk <- x[[7]]
           } else {
             nrowss <- dim(x[[7]])[1]
             alfas <- x[[7]] 
             pk <- x[[9]] 
           }
         }
         )

  rlabels <- rownames(alfas) 
  K <- ncol(alfas)  

  #### Obtendo as coordenadas por meio da análise de correspondência.
  N <- alfas
  P <- N/sum(N) 

  Umr <- matrix(1,nrow=ncol(P))
  Umc <- matrix(1,nrow=nrow(P)) 

  r <- as.vector(P%*%Umr)
  c <- as.vector(t(P)%*%Umc)

  S <- diag(1/sqrt(r))%*%(P-r%*%t(c))%*%diag(1/sqrt(c))

  svdd <- svd(S)

  aux_inertia <- svdd$d^2
  inertia <- aux_inertia[1:(length(aux_inertia)-1)]
  names(inertia) <- paste(round(inertia/sum(inertia)*100,2),
                          '%',
                          sep='')

  aux_rowcoordi <- diag(1/sqrt(r))%*%svdd$u 
  rowcoordi <- as.matrix(aux_rowcoordi[,-c(dim(aux_rowcoordi)[2])])
  colnames(rowcoordi) <- paste('Dim',
                               1:dim(rowcoordi)[2],
                               paste('(',
                                     names(inertia),
                                     ')',
                                     sep=''),
                               sep=' ')
  rownames(rowcoordi) <- rownames(N)

  aux_colcoordi <- diag(1/sqrt(c))%*%svdd$v
  colcoordi <- as.matrix(aux_colcoordi[,-c(dim(aux_colcoordi)[2])])
  colnames(colcoordi) <- paste('Dim',
                               1:dim(colcoordi)[2],
                               paste('(',
                                     names(inertia),
                                     ')',
                                     sep=''),
                               sep=' ')
  rownames(colcoordi) <- colnames(N)

  ################ Fim do algorítimo

  if(is.null(labels.budget)) labels.budget <- colnames(N)

  if(is.null(cex.budget)) cex.budget <- 1.2 

  if(is.null(col.budget)) col.budget <- rgb(0,0,0,0.3)

  if(is.null(pch.budget)) pch.budget <- 20 

  if(is.null(pos.budget)) pos.budget <- 3

  if(is.null(labels.points)) labels.points <- as.character(1:nrowss) 

  if(is.null(pch.points)) pch.points <- 10 

  if(is.null(col.points)) col.points <- 1

  if(is.null(xlab)) xlab <- colnames(rowcoordi)[dim][1]

  if(is.null(ylab)) ylab <- colnames(rowcoordi)[dim][2]

  if(is.null(pos.points)) pos.points <- 3

  plot(x    = rowcoordi[,dim],
       y    = NULL,
       type = 'n',
       xlim = xlim,
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       ...) 

  points(x = rowcoordi[,dim],
         pch = pch.points,
         col = col.points)

  text(x = rowcoordi[,dim],
       y = NULL,
       labels = labels.points,
       col = col.points,
       pos = pos.points,
       ...) 

  points(x = colcoordi[,dim],
         pch = pch.budget,
         col = col.budget)

  text(x = colcoordi[,dim],
       y = NULL,
       labels = labels.budget,
       cex = cex.budget,
       col = col.budget,
       pos = pos.budget,
       ...)

  abline(v = 0,
         lty = 3)

  abline(h = 0,
         lty = 3)

  if(is.null(args.legend)){

    thelabels <- paste(as.character(1:nrowss),
                       ' ',
                       rlabels,
                       sep="")

    args.2Kl <- list(x        = 'topleft',
                     legend   = thelabels,
                     text.col = col.points)

  } else {

    thelabels <- paste(as.character(1:nrowss),
                       ' ',
                       rlabels,
                       sep="")

    args.2Kl <- list(x        = 'topleft',
                     legend   = thelabels,
                     text.col = col.points)

    args.2Kl[names(args.legend)] <- args.legend     

  }       

  do.call('legend',
          args.2Kl) 

  coordenates <- list(rowcoordi,
                      colcoordi)

  invisible(coordenates)     

}
