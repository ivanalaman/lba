\name{summary.lba}
\alias{summary.lba.ls}
\alias{summary.lba.ls.fe}
\alias{summary.lba.ls.logit}
\alias{summary.lba.mle}
\alias{summary.lba.mle.fe}
\alias{summary.lba.mle.logit}

\title{
  Summary Method for \code{lba} objects.
}

\description{
  Returns (and prints) a summary list for objects of class \code{lba}, \code{lba.ls.fe}, \code{lba.ls.logit}, \code{lba.mle}, \code{lba.mle.fe}, and \code{lba.mle.logit}.
}

\usage{
\method{summary}{lba.ls}(object, digits = 2L, \dots)

\method{summary}{lba.ls.fe}(object, digits = 2L, \dots) 

\method{summary}{lba.ls.logit}(object, digits = 2L, \dots)  

\method{summary}{lba.mle}(object, digits = 2L, \dots)   

\method{summary}{lba.mle.fe}(object, digits = 2L, \dots)    

\method{summary}{lba.mle.logit}(object, digits = 2L, \dots)     
}

\arguments{
  \item{object}{A given object of the class \code{lba}, \code{lba.ls.fe}, \code{lba.mle.fe}, \code{lba.ls.logit} and \code{lba.mle.logit}.}
  \item{digits}{Number of decimal digits in the results.The default is 2.}      
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
summary(ex1)
 
}

\keyword{package}
\keyword{lba}
