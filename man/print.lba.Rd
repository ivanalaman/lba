\name{print.lba}
\alias{print.lba.ls}
\alias{print.lba.ls.fe}
\alias{print.lba.ls.logit}
\alias{print.lba.mle}
\alias{print.lba.mle.fe}
\alias{print.lba.mle.logit}

\title{
  Print Method for \code{lba} objects.
}

\description{
  Returns (and prints) a summary list for objects of class \code{lba.ls}, \code{lba.ls.fe}, \code{lba.ls.logit}, \code{lba.mle}, \code{lba.mle.fe}, and \code{lba.mle.logit}.
}

\usage{
\method{print}{lba.ls}(x, digits = 3L, ...)

\method{print}{lba.ls.fe}(x, digits = 3L, ...) 

\method{print}{lba.ls.logit}(x, digits = 3L, ...)  

\method{print}{lba.mle}(x, digits = 3L, ...)   

\method{print}{lba.mle.fe}(x, digits = 3L, ...)    

\method{print}{lba.mle.logit}(x, digits = 3L, ...)     
}

\arguments{
  \item{x}{A given object of the class \code{lba}, \code{lba.ls.fe}, \code{lba.mle.fe}, \code{lba.ls.logit} and \code{lba.mle.logit}.}
  \item{digits}{Number of decimal digits in the results. The default is 3.}      
  \item{\dots }{Further arguments (require by generic).}
}

\author{
  Enio G. Jelihovschi (\email{eniojelihovs@gmail.com})\cr
  Ivan Bezerra Allaman (\email{ivanalaman@gmail.com})\cr
}

\seealso{\code{\link{lba}}
}

\examples{
data('votB') 

# Using LS method (default) without constraint
# K = 2
ex1 <- lba(city ~ parties,
           votB,
           K = 2)
ex1
 
}

\keyword{package}
\keyword{lba}
