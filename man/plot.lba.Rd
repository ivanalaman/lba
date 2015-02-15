\name{plot.lba}
\alias{plot.lba}

\title{
 Plot lba objects
}

\description{
 S3 methods for \code{lba} objects.
 It's possible two types of visualisation: the lba type, suggested by van der Ark (1999) and correspondence analysis suggested by Jelihovschi (2011). 
}

\usage{
\method{plot}{lba}(x,
     budget.prop = TRUE,
     col.points  = NULL,
     col.lines   = NULL,# only to K = 2
     col.budget  = NULL,
     pch.points  = NULL,
     pch.budget  = NULL,
     lty.lines   = NULL,# only to K = 2
     lty.budget  = NULL,
     lwd.lines   = NULL,# only to K = 2
     lwd.budget  = NULL,
     legend      = TRUE,
     with.ml     = c('mix','lat'),
     type        =c("lba","corr"),
     ...) 
}

\arguments{
 \item{x}{A object of \code{lba} class.}
 \item{budget.prop}{A line representing the budget proportion. The default is \code{TRUE}.}
 \item{col.points}{A \code{vector} of colours representing the mixing parameters and latent components. The default is \code{NULL}.}
 \item{col.lines}{A \code{vector} of colours representing the mixing parameters and latent components. Only to K = 2. The default is \code{NULL}.}
 \item{col.budget}{A \code{vector} of colour representing the budget proportion. The default is \code{NULL}.}
 \item{pch.points }{A \code{vector} of plotting characters or symbols representing the mixing parameters and latent components. The default is \code{NULL}.}
 \item{pch.budget}{A \code{vector} of plotting characters or symbols representing the budget proportion. The default is \code{NULL}.}
 \item{lty.lines}{A \code{vector} of line types representing the mixing parameters and latent components. Only to K = 2. The default is \code{NULL}.}
 \item{lty.budget}{A \code{vector} of line types representing the budget proportion. The default is \code{NULL}.}
 \item{lwd.lines}{A \code{vector} of line width representing the mixing parameters and latent components. The default is \code{NULL}.}
 \item{lwd.budget}{A \code{vector} of line width representing the budget proportion. The default is \code{NULL}.}
 \item{legend}{A logical indicating whether a legend should be included. The default is \code{NULL}.}
 \item{with.ml}{What's parameters do you like to plot? The default is mixing parameters (\code{'mix'}).}
 \item{type}{The type of options graphical. The default is \code{'lba'}. See details.}
 \item{\dots}{Optional plotting parameters.}
}

\details{
 Plot.lba plots two types of graphics. The first, \code{(type = 'lba')}, is the one suggested at de Leeuw et all (1990) and at van der Ark (1999) thesis. In this type, it is only possible the graphical views for K = 2 and K = 3. When K = 2, the heads of the latent budgets are be connected by a (one dimensional) line segment. When K = 3, the heads of the latent budgets are the vertices of a triangle, and the plot is made with help of \code{triax.plot} function of \pkg{plotrix} package. In the second type, \code{plot.ca}, suggested at Jelihovschi et all (2011), the graphical display is made by making use of correspondence analysis graphics of the mixing parameters and latent components matrices. This is done with the function \code{ca} of \pkg{ca} package. In this case, a graphic display is possible for K >= 3.
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
 \code{\link[plotrix]{triax.plot}}, \code{\link[plotrix]{triax.points}}, \code{\link[ca]{plot.ca}}.
}

\examples{
data('votB')
ex1 <- lba(city ~ parties,
           data=votB,
           K = 2) 
plot(ex1)

# It's very simple. with colors!
plot(ex1,
     col.points=1:6,
     col.lines=1:6)

# I don't want your output. I want to provide my values and to use only your graphics. Ok!!!!
# Example: pag. 34, inner extreme solution
k1 <- c(.16,.68,.06,.45,.54,.30,.00,.45,1  ,.19,.62,.68,.57,.50,.30,
        .88,.46,.09,.62,.26,.53,.47,.58,.09,.26,.40,.70,.16,.71,.32)
k2 <- c(.84,.32,.94,.55,.46,.7 ,  1,.55,.00,.81,.38,.32,.43,.50,.7,
        .12,.54,.91,.38,.74,.47,.53,.42,.91,.74,.60,.30,.84,.29,.68)
my_mixing <- cbind(k1,k2)
rownames(my_mixing) <- paste(rep(c('H','C','R','Q','T','N','Se','Sk','L','V'),
                                 rep(3,10)),
                             c('+','=','-'),
                             sep = '')

my_budget <- matrix(c(.43,.57),
                    ncol=2)

colnames(my_budget) <- paste('LB',1:2,sep='')

my_data <- list(NA,
                NA,
                NA,
                NA,
                NA,
                my_mixing,
                NA,
                NA,
                my_budget)
class(my_data) = 'lba'

# I didn't like. I want equal of book. Pag. 41.

l_points <- c(2,24,17,
              rep(3,3),
              6,6,25,
              rep(4,3),
              1,16,21,
              rep(11,3),
              5,5,23,
              rep(8,3),
              rep(9,3),
              rep(13,3))

my_colors <- rep(c('red','blue','black'),
                 10)

plot(my_data,
     pch.points = l_points,
     col.points = my_colors,
     lty.budget = 2,
     pch.budget = 7,
     legend = FALSE)
legend(-0.1,
       1.2,
       rownames(my_mixing),
       pch = l_points,
       col = my_colors,
       xpd = TRUE,
       cex = 0.8,
       ncol = 10) # Is beautiful!! 
}

\keyword{lba}
\keyword{correpondence analysis}
\keyword{plot}