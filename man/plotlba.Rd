\name{plotlba}
\alias{plotlba}
\alias{plotlba.lba.1d}
\alias{plotlba.lba.2d}

\title{
  Plotlba objects using the approach suggested by van der Ark (1999).   
}

\description{
  S3 methods for \code{lba} objects. 
}

\usage{

\method{plotlba}{lba.1d}(x,
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

\method{plotlba}{lba.2d}(x,
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
}

\arguments{
  \item{x}{A object of \code{lba} class.}
  \item{height.line}{Is a vector with the lines height when \code{K = 2}.}
  \item{xlab}{A title for the x axis.}
  \item{ylab}{A title for the y axis.}
  \item{ylim}{The y limits of the plot.}
  \item{args.legend}{List of additional arguments to be passed to \code{legend}; names of the list are used as argument names. The default is \code{NULL}.}
  \item{axis.labels}{Labels for the three axes in the order left, right,
          bottom. Defaults to the column names.}
  \item{labels.points}{A character vector or expression specifying the text to be written. The default is \code{NULL}.}
  \item{col.points}{A \code{vector} of colour representing the points of the mixing parameters. The default is \code{par('fg')}.}
  \item{col.lines}{A \code{vector} of colour representing the lines of the mixing parameters. The default is \code{par('fg')}.}
  \item{lty.lines}{A \code{vector} of line types representing the mixing parameters. The default is \code{par('lty')}.} 
  \item{lwd.lines}{A \code{vector} of line width representing the mixing parameters. The default is \code{par('lwd')}.} 
  \item{pch.budget}{A \code{vector} of plotting characters or symbols representing the budget proportion. The default is \code{par('pch')}.} 
  \item{col.budget}{A \code{vector} of colour representing the budget proportion. The default is \code{par('fg')}.}
  \item{lty.budget}{A \code{vector} of line types representing the budget proportion. The default is \code{par('lty')}.}
  \item{lwd.budget}{A \code{vector} of line width representing the budget proportion. The default is \code{par('lwd')}.}
  \item{colline.budget}{The colors for line budget. The default is \code{par('fg')}.}
  \item{with.ml}{What's parameters do you like to plot? The default is mixing parameters (\code{'mix'}).}
  \item{\dots}{Other graphical parameters may also be passed as arguments to these functions.}
}

\details{
  The \code{plotlba} function, was suggested at de Leeuw et all (1990) and at van der Ark (1999) thesis. Those types of plots have only graphical views for \code{K = 2} and \code{K = 3}. When \code{K = 2}, either the latent budgets or the mixing parameters are displayed on a (one dimensional) line segment. When \code{K = 3}, either the latent budgets or the mixing parameters are displayed in a equilateral triangle using a barycentric coordinate system  where the budgets are represented by the vertices and the plot is made with help of \code{triax.plot} and \code{triax.points} function of \pkg{plotrix} package.
}

\references{
  de Leeuw, J., van der Heijden, P.G.M., and Verboon, P. 1990. "A latent time budget model". \emph{Statistica Neerlandica}. 44, 1, 1-21.

  van der Ark, A. L. 1999. \emph{Contributions to Latent Budget Analysis, a tool for the analysis of comositional data}. Ph.D. Thesis University of Utrecht.
}

\author{
  Enio G. Jelihovschi (\email{eniojelihovs@gmail.com})\cr
  Ivan Bezerra Allaman (\email{ivanalaman@gmail.com})
}

\seealso{
  \code{\link[plotrix]{triax.plot}}, \code{\link[plotrix]{triax.points}}.
}

\examples{
data('votB')

#K = 2
ex1 <- lba(parties ~ city,
           data=votB,
           K = 2) 

plotlba(ex1)

#It's very simple. With colors!
plotlba(ex1,
        col.points = 1:6,
        col.lines  = 1:6)        

#Add title in plot!
plotlba(ex1,
        main='Mixing parameters')

#Change budget proportion!
plotlba(ex1,
        pch.budget = 23,
        col.budget = 9,
        colline.budget = 8,
        lwd.budget = 2,
        lty.budget = 2)

#A little more!
plotlba(ex1,
        xlab = 'Lb2 -> Lb1',
        height.line = rep(0.5,6),
        lty.lines = 2,
        args.legend = list(ncol=3))


\dontrun{
#K = 3
data(MANHATAN)

tbm <- xtabs(value ~ socecon+health, 
             MANHATAN) 

ex2 <- lba(tbm,
           K = 3) 

plotlba(ex2)

#A little more!
plotlba(ex2,
        labels.points = rownames(tbm),
        col.points    = 2:7,
        args.legend   = list(plot=F))

plotlba(ex2,
        col.points     = 3:8,
        col.budget     = 2,
        pch.budget     = 20,
        lty.budget     = 2,
        lwd.budget     = 3,
        colline.budget = 3,
        axis.labels    = c('Lba1','Lba2','Lba3'))

}
}
\keyword{lba}
\keyword{plot}
