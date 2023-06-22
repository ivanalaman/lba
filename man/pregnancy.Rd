\name{pregnancy}
\alias{pregnancy}
\docType{data}

\title{Pregnancy-Related Mortality in California}

\description{
The \code{pregnancy} matrix has 16 rows and 5 columns. The raw data refers to California pregnancy-related deaths from 2002-2005.
}

\format{
 This matrix contains the following columns:
 \describe{
  \item{Pre.E}{Preeclampsia/eclampsia}
  \item{OH}{Obstetric hemorrhage}
  \item{CVD}{Cardiovascular diseases}
  \item{DVTPE}{Deep vein thrombosis - pulmonary embolism}
  \item{AFE}{Amniotic fluid embolism}
 }
 The rows refers to:
 \describe{
  \item{Hfob}{Hispanic,foreign-born}
  \item{Husb}{Hispanic, us-born}
  \item{Wnh}{White, non-hispanic}
  \item{Bnh}{Black, non-hispanic}
  \item{$<30$b}{Maternal age}     
  \item{$30-40$b}{Maternal age}
  \item{$>40$b}{Maternal age}
  \item{$1$}{Parity}
  \item{$2-4$}{Parity}
  \item{$5+$}{Parity}
  \item{$<30$a}{Maternal age}
  \item{$30-40$a}{Maternal age}
  \item{$>40$a}{Maternal age}
  \item{$<32$w}{Gestational age at delivery}
  \item{$32-36$w}{Gestational age at delivery}
  \item{$>37$w}{Gestational age at delivery}
 }
}

\source{
Main, E. K.; et al. Pregnancy-Related Mortality in California: Causes, Characteristics, and Improvement Opportunities. \emph{OBSTETRICS and GYNECOLOGY}. \bold{125}, 938--947.
}

\references{
Main, E. K.; et al. Pregnancy-Related Mortality in California: Causes, Characteristics, and Improvement Opportunities. \emph{OBSTETRICS and GYNECOLOGY}. \bold{125}, 938--947.
}

\usage{
data(pregnancy)
}

\keyword{dataset}
