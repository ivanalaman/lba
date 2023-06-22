plotcorr.lba.3d <- function(x,
                            rgl.use        = FALSE,
			    dim            = c(1,2,3), #only K >= 3
                            xlim           = NULL,
                            ylim           = NULL,
                            zlim           = NULL,
			    xlab           = NULL,
			    ylab           = NULL,
			    zlab           = NULL,
			    args.legend    = NULL, #only rgl.use=FALSE
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

  if(is.null(labels.budget)) labels.budget <- colnames(N)#ok

  if(is.null(cex.budget)) cex.budget <- 1.2#ok 

  if(is.null(col.budget)) col.budget <- rgb(0,0,0,0.3)#ok

  if(is.null(pch.budget)) pch.budget <- 20#ok 

  if(is.null(pos.budget)) pos.budget <- 3#ok

  if(is.null(labels.points)) labels.points <- as.character(1:nrowss)#ok 

  if(is.null(pch.points)) pch.points <- 2#ok 

  if(is.null(pos.points)) pos.points <- 3

  if(is.null(col.points)) col.points <- 1

  if(is.null(xlab)) xlab <- colnames(rowcoordi)[dim][1]

  if(is.null(ylab)) ylab <- colnames(rowcoordi)[dim][2] 

  if(is.null(zlab)) zlab <- colnames(rowcoordi)[dim][3]  

  if(!rgl.use){

    rowcoordis <- rowcoordi[,dim]
    colcoordis <- colcoordi[,dim]

  gg <- scatterplot3d(rowcoordis,
                      xlim  = xlim,
                      ylim  = ylim,
                      zlim  = zlim,
		      xlab  = xlab,
		      ylab  = ylab,
		      zlab  = zlab,
                      pch   = pch.points, 
                      color = col.points,
		      ...)

  gg$points3d(colcoordis,
	      pch = pch.budget,
	      ...)
  
  gg$points3d(seq(min(rowcoordis[,1]),
                  max(rowcoordis[,1]),
                  l = dim(rowcoordis)[1]),
              rep(0,dim(rowcoordis)[1]),
              rep(0,dim(rowcoordis)[1]),
              type = 'l',
              lty = 3)

  gg$points3d(rep(0,dim(rowcoordis)[1]),
              seq(min(rowcoordis[,2]),
                  max(rowcoordis[,2]),
                  l = dim(rowcoordis)[1]), 
              rep(0,dim(rowcoordis)[1]),
              type = 'l',
              lty = 3)

  gg$points3d(rep(0,dim(rowcoordis)[1]),
              rep(0,dim(rowcoordis)[1]),
              seq(min(rowcoordis[,3]),
                  max(rowcoordis[,3]),
                  l = dim(rowcoordis)[1]),  
              type = 'l',
              lty = 3)

  row.coords <- gg$xyz.convert(rowcoordi[,dim])
  col.coords <- gg$xyz.convert(colcoordi[,dim]) 

  text(x = row.coords$x, 
       y = row.coords$y, 
       labels = labels.points,
       pos = pos.points,
       color = col.points,
       offset = 0.5,
       ...)

  text(x      = col.coords$x, 
       y      = col.coords$y, 
       labels = labels.budget,
       cex    = cex.budget,
       color    = col.budget,
       pos    = pos.budget,
       offset = 0.5)

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

  } else {
  
    plot3d(rowcoordi[,dim],
           xlab  = xlab,
           ylab  = ylab,
           zlab  = zlab,
           pch   = pch.points, 
           color   = col.points,
           ...)

    aux1_rowcoordi <- rowcoordi[,dim] 
    text3d(x = aux1_rowcoordi[,1],
           y = aux1_rowcoordi[,2],
           z = aux1_rowcoordi[,3]+0.3,
           texts = labels.points,
           color = col.points,
           ...)

    pch3d(colcoordi[,dim],
          pch = pch.budget,
          color = col.budget,
          ...)

    aux1_colcoordi <- colcoordi[,dim]
    text3d(x = aux1_colcoordi[,1],
           y = aux1_colcoordi[,2],
           z = aux1_colcoordi[,3]+0.3,
           texts = labels.budget,
           color = col.budget,
           ...)


  }

  coordenates <- list(rowcoordi[,dim],
                      colcoordi[,dim])

  invisible(coordenates)      
}
