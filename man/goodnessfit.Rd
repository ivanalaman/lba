\name{goodnessfit}
\alias{goodnessfit}

\title{Goodness of Fit results for Latent Budget Analysis}

\description{The goodness of fit results assesses how well the model fits the data. It consists of measures of the resemblance between the observed and the expected data, and the parsimony of the model.
}

\usage{
goodnessfit(object)
}

\arguments{
 \item{object}{An object of one of following classes: \code{lba.ls}, \code{lba.ls.fe}, \code{lba.ls.logit}, \code{lba.mle}, \code{lba.mle.fe}, \code{lba.mle.logit}}
}

\value{
  The \code{goodnessfit} function to the method \code{lba.mle}, \code{lba.mle.fe} and \code{lba.mle.logit} returns a list with the slots:
      
 \item{dfdb     }{Degrees of freedom of the base model}
 \item{dfd      }{Degrees of freedom of the full model}
 \item{G2b      }{Likelihood ratio statistic of the base model}
 \item{G2       }{Likelihood ratio statistic of the full model}
 \item{chi2b    }{Chi-square statistic of the base model}
 \item{chi2     }{Chi-square statistic of the full model}
 \item{proG1    }{P-value of likelihood ratio statistic of the base model}
 \item{proG     }{P-value of likelihood ratio statistic of the full model}
 \item{prochi1  }{P-value of chi-square statistic of the base model}
 \item{prochi   }{P-value of chi-square statistic of the full model}
 \item{AICb     }{AIC criteria of the base model}
 \item{AICC     }{AIC criteria of the full model}
 \item{BICb     }{BIC criteria of the base model}
 \item{BICC     }{BIC criteria of the full model}
 \item{CAICb    }{CAIC criteria of the base model}
 \item{CAIC     }{CAIC criteria of the full model}
 \item{delta1   }{Normed fit index}
 \item{delta2   }{Normed fit index modified}
 \item{rho1     }{Bollen index}
 \item{rho2     }{Tucker-Lewis index}
 \item{RSS1     }{Residual sum of square of the base model}
 \item{RSS      }{Residual sum of square of the full model}
 \item{impRSS   }{Improvement of RSS}
 \item{impPB    }{Improvement per budget}
 \item{impDF    }{Average improvement per degree of freedom}
 \item{D1       }{Index of dissimilarity of the base model}
 \item{D        }{Index of dissimilarity of the full model}
 \item{pccb     }{Proportion of correctly classified data of the base model}
 \item{pcc      }{Proportion of correctly classified data of the full model}
 \item{impD     }{Improvement of proportion of correctly classified data}
 \item{impPCCB  }{Improvement of Proportion of correctly classified data per budget}
 \item{AimpPCCDF}{Average improvement of Proportion of correctly classified data per degree of freedom}
 \item{mad1     }{Mean angular deviation of the base model}
 \item{madk     }{Mean angular deviation of the full model}
 \item{impMad   }{Improvement mean angular deviation}
 \item{impPBsat }{Improvement mean angular deviation per budget}
 \item{impDFsat }{Average improvement mean angular deviation per degree of freedom}

 The \code{goodnessfit} function to the method \code{lba.ls}, \code{lba.ls.fe} and \code{lba.ls.logit} returns a list with the slots:

 \item{dfdb     }{Degrees of freedom of the base model}
 \item{dfd      }{Degrees of freedom of the full model}
 \item{RSS1     }{Residual sum of square of the base model}
 \item{RSS      }{Residual sum of square of the full model}
 \item{impRSS   }{Improvement of RSS}
 \item{impPB    }{Improvement per budget}
 \item{impDF    }{Average improvement per degree of freedom}
 \item{wRSS1    }{Weighted residual sum of square of the base model}
 \item{wRSS     }{Weighted residual sum of square of the full model}
 \item{impwRSS  }{Improvement of wRSS}
 \item{D1       }{Index of dissimilarity of the base model}
 \item{D        }{Index of dissimilarity of the full model}
 \item{pccb     }{Proportion of correctly classified data of the base model}
 \item{pcc      }{Proportion of correctly classified data of the full model}
 \item{impD     }{Improvement of proportion of correctly classified data}
 \item{impPCCB  }{Improvement of Proportion of correctly classified data per budget}
 \item{AimpPCCDF}{Average improvement of Proportion of correctly classified data per degree of freedom}
 \item{mad1     }{Mean angular deviation of the base model}
 \item{madk     }{Mean angular deviation of the full model}
 \item{impMad   }{Improvement mean angular deviation}
 \item{impPBsat }{Improvement mean angular deviation per budget}
 \item{impDFsat }{Average improvement mean angular deviation per degree of freedom}
}

\references{
   Agresti, Alan. 2002. \emph{Categorical Data Analysis, second edition}. Hoboken: John Wiley \& Sons.

   van der Ark, A. L. 1999. \emph{Contributions to Latent Budget Analysis, a tool for the analysis of compositional data}. Ph.D. Thesis University of Utrecht.
}

\note{For a detailed and complete discussion about goodness of fit results for latent budget analysis, see van der Ark 1999.
}

\seealso{
 \code{\link{print.goodnessfit}}, \code{\link{lba}}
}

\examples{
data('votB')

# Using LS method (default) without constraint
# K = 2
ex1 <- lba(city ~ parties,
           votB,
           K = 2)

gx1 <- goodnessfit(ex1)
gx1

# Using MLE method without constraint
# K = 2
exm <- lba(city ~ parties,
           votB,
           K = 2,
           method='mle')

gxm <- goodnessfit(exm)
gxm
}