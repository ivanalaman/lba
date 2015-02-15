\name{votB}
\alias{votB}
\docType{data}

\title{Voting Behaviour in Netherlands}

\description{

The \code{votB} data frame has 8971 rows and 2 columns.
The raw data refers to the type of the city and the political party which each 
participant voted for in the 1986 general elections in the Netherlands.
}

\format{
 This data frame contains the following columns:
 \describe{
  \item{city}{
   A factor with levels:
   \code{co} Commuter;
   \code{lx} Large city;
   \code{mc} Middle large city;
   \code{ri} Rural industrialised;
   \code{ru} Rural;
   \code{sc} Small city.
  }
  \item{parties}{
   A factor with levels:
   \code{cda  } Christian democrats;
   \code{d66  } Democrats;
   \code{left } Other left-wing parties;
   \code{Pvda } Labor party;
   \code{right} Other right-wing parties;
   \code{vvd  } Liberals.
  }
 }
}

\source{
 Statistics Netherlands (1987). Statistiek der verkiezingen 1986. Tweede Kamer der Staten-Generaal 21 mei 1996. [Statistics of the elections of the Lower House, May 21-th 1996]. The Hague: Staatsuitgeverij.
}

\references{
 van der Ark, A. L. 1999. \emph{Contributions to Latent Budget Analysis, a tool for the analysis of compositional data}. Ph.D. Thesis University of Utrecht. 
}

\usage{
votB
}

\keyword{dataset}