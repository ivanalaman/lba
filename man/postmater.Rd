\name{postmater}
\alias{postmater}
\docType{data}

\title{Environmental problems and Cancer}

\description{

The \code{postmater} data frame has 13 rows and 8 columns. The raw data refers to a political and social survey across Europe that is conducted twice a year.
}

\format{
 This data frame contains the following columns:
 \describe{
  \item{country}{
   A factor with levels:
   \code{F} France;
   \code{B} Belgium;
   \code{NL} Netherlands;
   \code{D} Germany;
   \code{I} Italy;
   \code{L} Luxembourg.
   \code{DK} Denmark. 
   \code{IRL} Ireland. 
   \code{GB} Great Britain. 
   \code{NIRL} Northern Ireland. 
   \code{GR} Greece. 
   \code{E} Spain. 
   \code{P} Portugal.  
  }
  \item{m..}{
    The absolute frequencies of materialist factor in the respect country. The degree of ranking of this index is ++.
  }
  \item{m.}{
    The absolute frequencies of materialist factor in the respect country. The degree of ranking of this index is +.
  } 
  \item{m}{
    The absolute frequencies of materialist factor in the respect country. The degree of ranking of this index is below of the "m.".
  } 
  \item{m_pm}{
    The absolute frequencies of materialist/post-materialist factor in the respect country. The degree of ranking of this index is below of the "m".
  }
  \item{pm}{
    The absolute frequencies of post-materialist factor in the respect country. The degree of ranking of this index is below of the "m_pm".
  } 
  \item{pm.}{
    The absolute frequencies of post-materialist factor in the respect country. The degree of ranking of this index is below of the "pm".
  }  
  \item{pm..}{
    The absolute frequencies of post-materialist factor in the respect country. The degree of ranking of this index is below of the "pm.".
  }  
 }
}

\source{
Reif, K., and Melich, A. (1990). Euro-Barometer 29: Environmental problems and Cancer, March-April 1988. Ann Arbor: Inter-university Consortium for Political and Social Research.
}

\references{
 van der Ark, A. L. 1999. \emph{Contributions to Latent Budget Analysis, a tool for the analysis of compositional data}. Ph.D. Thesis University of Utrecht. 
}

\usage{
postmater
}

\keyword{dataset} 
