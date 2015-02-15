\name{print.goodnessfit}
\alias{print.goodnessfit}

\title{
  Print Method for goodnessfit objects.
}

\description{
  Returns (and prints) a summary list for goodnessfit objects.
}

\usage{
\method{print}{goodnessfit}(x, digits = 2L, \dots)
}

\arguments{
  \item{x}{A given object of the class \code{goodnessfit.lba.ls} and \code{goodnessfit.lba.mle}.}
  \item{digits}{Number of decimal digits in the results.}      
  \item{\dots }{Potential further arguments (require by generic).}
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
exm 
}

\keyword{package}
\keyword{goodnessfit}