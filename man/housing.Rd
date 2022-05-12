\name{housing}
\alias{housing}
\docType{data}

\title{The Satisfaction with Housing Conditions Study}

\description{
The \code{housing} data frame has 72 rows and 5 columns.
The observations were obtained from an investigation of Satisfaction with housing conditions carried out by the Danish Building Research Institute and the Danish Institute of Mental Health Research.
}

\format{
 This data frame contains the following columns:
 \describe{
  \item{Housing}{
   A factor with levels:
   \code{Apartment};
   \code{Atrium};
   \code{Terraced};
   \code{Tower}.
  }
  \item{Influence}{
   A factor with levels:
   \code{hi};
   \code{low};
   \code{med}.
  }
  \item{Contact}{
   A factor with levels:
   \code{high};
   \code{low}.
  }
 \item{Satisfaction}{
   A factor with levels:
   \code{hi};
   \code{low};
   \code{med}.
  } 
  \item{value}{
   The absolute frequencies of which factor. 
  }
 }
}

\source{
 Madsen, M. (1976) Statistical analysis of multiple contingency tables: Two examples. \emph{Scandinavian Journal of Statistics} \bold{3}, 97--106.
}

\references{
 van der Ark, A. L. 1999. \emph{Contributions to Latent Budget Analysis, a tool for the analysis of compositional data}. Ph.D. Thesis University of Utrecht. 
}

\usage{
data(housing)
}

\keyword{datasets}
