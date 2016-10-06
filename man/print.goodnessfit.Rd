\name{print.goodnessfit}
\alias{print.goodnessfit.lba.ls}
\alias{print.goodnessfit.lba.ls.fe}
\alias{print.goodnessfit.lba.ls.logit}
\alias{print.goodnessfit.lba.mle}
\alias{print.goodnessfit.lba.mle.fe}
\alias{print.goodnessfit.lba.mle.logit}

\title{
  Print Method for goodnessfit objects.
}

\description{
  Returns (and prints) a summary list for goodnessfit objects.
}

\usage{
\method{print}{goodnessfit.lba.ls}(x, digits=3L, \dots)

\method{print}{goodnessfit.lba.ls.fe}(x, digits=3L, \dots) 

\method{print}{goodnessfit.lba.ls.logit}(x, digits=3L, \dots)  

\method{print}{goodnessfit.lba.mle}(x, digits=3L, \dots)   

\method{print}{goodnessfit.lba.mle.fe}(x, digits=3L, \dots)    

\method{print}{goodnessfit.lba.mle.logit}(x, digits=3L, \dots)     

}

\arguments{
  \item{x}{A given object of the class \code{goodnessfit.lba.ls}, \code{goodnessfit.lba.ls.fe}, \code{goodnessfit.lba.ls.logit}, \code{goodnessfit.lba.mle.fe}, \code{goodnessfit.lba.mle.logit} and \code{goodnessfit.lba.mle}.}
  \item{digits}{A non-null value for digits specifies the minimum number of significant digits to be printed in values. The default is 3.}
  \item{\dots }{Further arguments (require by generic).}
}

\author{
  Enio G. Jelihovschi (\email{eniojelihovs@gmail.com})\cr
  Ivan Bezerra Allaman (\email{ivanalaman@gmail.com})\cr
}

\seealso{
  \code{\link{summary.goodnessfit.lba.ls}}, \code{\link{summary.goodnessfit.lba.mle}}
  }

\examples{
data('votB')

# Using LS method (default) without constraint
# K = 2
ex1 <- lba(city ~ parties,
           votB,
           K = 2)
exm <- goodnessfit(ex1)
exm

}

\keyword{package}
\keyword{goodnessfit}
