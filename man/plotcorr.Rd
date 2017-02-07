\name{plotcorr}
\alias{plotcorr}
\alias{plotcorr.lba.1d}
\alias{plotcorr.lba.2d} 
\alias{plotcorr.lba.3d}  

\title{
  Plot lba objects using the correspondence analysis approach as suggested by Jelihovschi (2011).
}

\description{
  S3 methods for \code{lba} objects. 
}

\usage{
plotcorr(x, ...)

\method{plotcorr}{lba.1d}(x,
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

\method{plotcorr}{lba.2d}(x,
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

\method{plotcorr}{lba.3d}(x,
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

}

\arguments{
  \item{x}{A object of \code{lba} class.}
  \item{dim}{The dimention to be plotted. The default is \code{c(1,2)} to \code{K = 2} and \code{c(1,2,3)} to \code{K = 3}.}
  \item{xlim}{The x limits (x1, x2) of the plot.}
  \item{ylim}{The y limits of the plot.}
  \item{zlim}{The z limits of the plot.}
  \item{xlab}{A label for the x axis, defaults to a description of "x".} 
  \item{ylab}{A label for the y axis, defaults to a description of "y".} 
  \item{zlab}{A label for the z axix, defaults to a description of "z".} 
  \item{rgl.use}{A logical value. If \code{TRUE} the 3d scatter will be done with the \code{rgl} environment, in another way the \code{scatterplot3d} will be used.}
  \item{metrics}{Logical. If \code{TRUE} (default), the \code{radius} is plotted.}
  \item{radius}{A arbitrary number to choose the groups. The default is \code{0.5}. See details.}
  \item{col.points}{The color points to be used, possibly vectors. The default is \code{NULL}. See datails.} 
  \item{height.points}{Points label height in relation to the y-coordinate.The default is \code{NULL}.} 
  \item{labels.points}{A character vector or expression specifying the _text_ to be written. The default is \code{NULL}.}
  \item{pch.points}{A symbols to use. O default is \code{NULL}.}
  \item{pos.points}{A position specifier for the text.  If specified this overrides any "adj" value given.  Values of "1", "2", "3" and "4", respectively indicate positions below, to the left of, above and to the right of the specified coordinates.}
  \item{args.legend}{List of additional arguments to be passed to \code{legend}; names of the list are used as argument names. Only used if \code{K=2}. The default is \code{NULL}.}
  \item{pch.budget}{A symbols to use. O default is \code{NULL}.}
  \item{pos.budget}{A position specifier for the text.  If specified this overrides any "adj" value given.  Values of "1", "2", "3" and "4", respectively indicate positions below, to the left of, above and to the right of the specified coordinates.}
  \item{height.budget}{Budget label height in relation to the y-coordinate.The default is \code{NULL}.}
  \item{labels.budget}{A character vector or expression specifying the _text_ to be written. The default is \code{NULL}.}
  \item{cex.budget}{The size of text. The default is \code{NULL}.}
  \item{col.budget}{The color budget to be used, possibly vectors. The default is \code{NULL}.}
  \item{with.ml}{Vector of two character strings specifying the parameters of the plot. Set \code{"mix"} to plot the mixing parameters and \code{"lat"} to plot the latent components. The default is \code{"mix"}.}
  \item{\dots}{Further graphical parameters.}
}

\details{
  The \code{plotcorr} suggested by Jelihovschi et all (2011), has a graphical display which uses the correspondence analysis graphics of the mixing parameters and latent components matrices. In this case, a graphic display is possible for \code{K >= 2}.

  The argument \code{radius} was featured in order to help the user as he or she needs do decide which are the points belonging to a certain latent budget. Only the points to the right or left of LB1 and LB2 but always towards the center of the graphic (the zero of x axis) were taken in account, since those in opposite direction automatically belong to the closest latent budget. this argument only works for \code{K = 2}. It's should be of size two. 

 The argument \code{col.points} takes in account the argument \code{radius} in order to color the groups which either belong or not to a certain budget, therefore, the size of the vector of this argument must be equal to the number of formed groups.

}

\references{
  de Leeuw, J., van der Heijden, P.G.M., and Verboon, P. 1990. "A latent time budget model". \emph{Statistica Neerlandica}. 44, 1, 1-21.

  Jelihovschi, E.G., Alves, R.R., and Correa, F.M. 2011. \emph{Interacting latent budget analysis and correspondence analysis to analyze beauty salon management data}. Biometric Brazilian Journal, 29, 657-673.

  van der Ark, A. L. 1999. \emph{Contributions to Latent Budget Analysis, a tool for the analysis of comositional data}. Ph.D. Thesis University of Utrecht.
}

\author{
  Enio G. Jelihovschi (\email{eniojelihovs@gmail.com})\cr
  Ivan Bezerra Allaman (\email{ivanalaman@gmail.com})
}

\seealso{
  \code{\link[graphics]{plot.default}}, \code{\link[scatterplot3d]{scatterplot3d}}, \code{\link[rgl]{plot3d}}.
}

\examples{
data('votB')
K = 2
ex1 <- lba(parties ~ city,
           data=votB,
           K = 2) 
plotcorr(ex1)

#It's very simple. with colors!
plotcorr(ex1,
         col.points = 3:5,
         col.budget = c(5,3))

#Changing radius!
plotcorr(ex1,
         radius = rep(0.7,2))

#Without metrics!
plotcorr(ex1,
         metrics = FALSE)

#Change legend options!
plotcorr(ex1,
         args.legend = list(ncol=3)) 

#Change height points!
plotcorr(ex1,
         height.points = rep(-0.1,6))


\dontrun{
#K = 3
K = 3
ex2 <- lba(parties ~ city,
           data=votB,
           K = 3)

plotcorr(ex2)

#Change budget options
plotcorr(ex2,
         pch.budget = 5,
         col.budget = 2,
         labels.budget = c('lba1','lba2','lba3'))

#Change points options
plotcorr(ex2,
         pch.points = 20,
         col.points = 4,
         labels.points = rownames(ex2$Aoi),
         args.legend = list(plot=FALSE))

#Coloring the groups
plotcorr(ex2,
         col.points = c(1,2,2,3,3,2),
         col.budget = c(3,1,2),
         args.legend = list(ncol=3))

#K = 4
K = 4
data(postmater)
new_post <- as.matrix(postmater[,-1])
row.names(new_post) <- postmater[,1]

ex3 <- lba(new_post,
           K = K)

plotcorr(ex3)



#A bit didatic!
plotcorr(ex3,
         args.legend = list(x = -2.5,
                            y = 5.5,
                            xpd=TRUE,
                            ncol=5))

#Dynamic? Yes, you can!
plotcorr(ex3,
     rgl.use = TRUE)
}
}

\keyword{lba}
\keyword{correpondence analysis}
\keyword{plot}

