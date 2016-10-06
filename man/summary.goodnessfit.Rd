\name{summary.goodnessfit}
\alias{summary.goodnessfit.lba.ls}
\alias{summary.goodnessfit.lba.ls.fe}
\alias{summary.goodnessfit.lba.ls.logit}
\alias{summary.goodnessfit.lba.mle}
\alias{summary.goodnessfit.lba.mle.fe}
\alias{summary.goodnessfit.lba.mle.logit}

\title{
  Summary Method for goodnessfit objects.
}

\description{
  Returns (and prints) a summary list for goodnessfit objects.
}

\usage{
\method{summary}{goodnessfit.lba.ls}(object, digits = 2L, \dots)

\method{summary}{goodnessfit.lba.ls.fe}(object, digits = 2L, \dots)
 

\method{summary}{goodnessfit.lba.ls.logit}(object, digits = 2L, \dots)
  

\method{summary}{goodnessfit.lba.mle}(object, digits = 2L, \dots)
   

\method{summary}{goodnessfit.lba.mle.fe}(object, digits = 2L, \dots)
    

\method{summary}{goodnessfit.lba.mle.logit}(object, digits = 2L, \dots)     
}

\arguments{
  \item{object}{A given object of the class \code{goodnessfit.lba.ls} and \code{goodnessfit.lba.mle}.}
  \item{digits}{Number of decimal digits in the results. The default is 2.}      
  \item{\dots }{Further arguments (require by generic).}
}

\author{
  Enio G. Jelihovschi (\email{eniojelihovs@gmail.com})\cr
  Ivan Bezerra Allaman (\email{ivanalaman@gmail.com})\cr
}

\seealso{\code{\link{goodnessfit}}
}

\examples{
data('votB')

# Using LS method (default) without constraint
# K = 2
ex1 <- lba(city ~ parties,
           votB,
           K = 2)
exm <- goodnessfit(ex1)
summary(exm)

}

\keyword{package}
\keyword{goodnessfit}
