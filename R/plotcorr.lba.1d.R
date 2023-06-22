plotcorr.lba.1d <- function(x,
                         xlim           = NULL,
                         ylim           = NULL,
                         xlab           = NULL,
                         ylab           = NULL,
                         metrics        = TRUE,
                         radius         = rep(0.5,2),
                         col.points     = NULL,
                         height.points  = NULL,     
                         labels.points  = NULL,
                         pch.points     = NULL,
                         pos.points     = NULL,
                         args.legend    = NULL, 
                         height.budget  = NULL,    
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

  #### Obtendo as coordenadas por meio da analise de correspondencia.
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

  ################ Fim do algoritimo

  if(is.null(labels.budget)) labels.budget <- colnames(N)

  if(is.null(cex.budget)) cex.budget <- 1.2 

  if(is.null(col.budget)) col.budget <- rgb(0,0,0,0.3)

  if(is.null(pch.budget)) pch.budget <- 20 

  if(is.null(pos.budget)) pos.budget <- 3

  if(is.null(labels.points)) labels.points <- as.character(1:nrowss) 

  if(is.null(pch.points)) pch.points <- 10 

  if(is.null(pos.points)) pos.points <- 3 

  if(colcoordi[1,1] > colcoordi[2,1]){

    colcoor <- -colcoordi
    rowcoor <- -rowcoordi

  } else {

    colcoor <- colcoordi
    rowcoor <- rowcoordi

  }

  ### Calculo das distancias

  groups_lb1 <- as.data.frame(apply(rowcoor,
                                    2,
                                    function(x) ifelse(x < colcoor[1,],
                                                       'lb1',
                                                       'nenhum')))
  groups_lb2 <- as.data.frame(apply(rowcoor,
                                    2,
                                    function(x) ifelse(x > colcoor[2,],
                                                       'lb2',
                                                       'nenhum')))

  dist_lb1 <- as.data.frame(apply(rowcoor, 
                                  2, 
                                  function(x) abs(x - colcoor[1,]))) 

  dist_lb2 <- as.data.frame(apply(rowcoor, 
                                  2, 
                                  function(x) abs(x - colcoor[2,]))) 

  if(length(radius) != 2) stop('The size should be two!') 
  in_group <- radius

  groups_lb11 <- apply(dist_lb1, 
                       2, 
                       function(x) ifelse(x < in_group[1], 
                                          'lb1',
                                          'nenhum')) 

  groups_lb22 <- apply(dist_lb2, 
                       2, 
                       function(x) ifelse(x < in_group[2], 
                                          'lb2',
                                          'nenhum')) 

  groups_aux <- cbind(groups_lb1,
                      groups_lb11,
                      groups_lb2,
                      groups_lb22)

  groups <- apply(groups_aux, 
                  1, 
                  function(x) ifelse(any(x=='nenhum'),
                                     x[x!='nenhum'],
                                     'nenhum'))
  groups[is.na(groups)] <- 'nenhum'

  n_groups <- unique(groups)

  if(!is.null(col.points)){
    if(length(col.points) != length(n_groups)) stop('The size of your color should be the same size your group!')

    col_aux <- col.points

  } else{

    col_aux <- 1:length(n_groups) 

  }

  mgsub <- function(pattern, replacement, x, ...) {
    if (length(pattern)!=length(replacement)) {
      stop("pattern and replacement do not have the same length.")
    }
    result <- x
    for (i in 1:length(pattern)) {
      result <- gsub(pattern[i], replacement[i], result, ...)
    }
    result
  }

  col_new <- as.numeric(mgsub(n_groups,
                              col_aux,
                              groups))

  #### Fim do algoritimo das distancias

  mia <- min(rowcoor, 
             colcoor)
  maa <- max(rowcoor, 
             colcoor)

  lin <- mia+(mia*0.3)
  lsu <- maa+(maa*0.3)

  if(is.null(height.points)) height.points <- seq(-0.5,0.5,l=nrowss) 

  if(is.null(xlim)) xlim <- seq(lin,lsu,l = 2)

  if(is.null(ylim)) ylim <- c(-1.5,1.5)

  if(is.null(xlab)) xlab <- ''

  if(is.null(ylab)) ylab <- ''

  if(is.null(height.budget)) height.budget <- rep(0,2)

  coor_points <- matrix(c(rowcoor,
                          height.points),
                        ncol = 2)

  plot(x    = coor_points,
       y    = NULL,
       type = 'n',
       axes = F,
       xlim = xlim,
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       ...)

  axis(1,
       ...)
  box()

  points(x = colcoor,
         y = height.budget,
         pch = pch.budget,
         col = col.budget,
         ...)

  text(x = colcoor,
       y = height.budget,
       labels = labels.budget,
       cex = cex.budget,
       col = col.budget,
       pos = pos.budget)

  points(x = coor_points,
         y = NULL,
         pch = pch.points,
         col = col_new,
         ...)

  text(x = coor_points,
       y = NULL,
       labels = labels.points,
       col = col_new,
       pos = pos.points,
       ...)

  if(metrics){ 

    segments(x0  = colcoor[1,],
             y0  = height.budget[1],
             x1  = colcoor[1,]+radius[1],
             y1  = height.budget[1],
             lty = 2
             )

    text(x = median(c(colcoor[1,],colcoor[1,]+radius[1])),
         y = height.budget[1]-0.1,
         labels = paste(radius[1],
                        'un.'),
         cex = 0.8)

    segments(x0  = colcoor[2,],
             y0  = height.budget[1],
             x1  = colcoor[2,]-radius[2],
             y1  = height.budget[1],
             lty = 2
             ) 

    text(x = median(c(colcoor[2,],colcoor[2,]-radius[2])),
         y = height.budget[1]-0.1,
         labels = paste(radius[2],
                        'un.'),
         cex = 0.8)   
  }

  if(is.null(args.legend)){

    thelabels <- paste(as.character(1:nrowss),
                       ' ',
                       rlabels,
                       sep="")

    args.2Kl <- list(x        = 'topleft',
                     legend   = thelabels,
                     text.col = col_new)

  } else {

    thelabels <- paste(as.character(1:nrowss),
                       ' ',
                       rlabels,
                       sep="")

    args.2Kl <- list(x        = 'topleft',
                     legend   = thelabels,
                     text.col = col_new)

    args.2Kl[names(args.legend)] <- args.legend     

  }       

  do.call('legend',
          args.2Kl) 

  coordenates <- list(coor_points,
                      colcoor)

  invisible(coordenates) 

} 
