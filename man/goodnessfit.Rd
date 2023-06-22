\name{goodnessfit}
\alias{goodnessfit}
\alias{goodnessfit.default}
\alias{goodnessfit.lba.ls}
\alias{goodnessfit.lba.ls.fe}
\alias{goodnessfit.lba.ls.logit}
\alias{goodnessfit.lba.mle}
\alias{goodnessfit.lba.mle.fe}
\alias{goodnessfit.lba.mle.logit}

\title{Goodness of Fit results for Latent Budget Analysis}

\description{The goodness of fit results assesses how well the model fits the data. It consists of measures of the resemblance between the observed and the expected data, and the parsimony of the model.
}

\usage{
goodnessfit(object,...)

## S3 methods
\method{goodnessfit}{default}(object, \dots)

\method{goodnessfit}{lba.ls}(object, \dots)
    
\method{goodnessfit}{lba.ls.fe}(object, \dots)    
    
\method{goodnessfit}{lba.ls.logit}(object, \dots)
    
\method{goodnessfit}{lba.mle}(object, \dots)    

\method{goodnessfit}{lba.mle.fe}(object, \dots)    
 
\method{goodnessfit}{lba.mle.logit}(object, \dots) 
}

\arguments{
 \item{object}{An object of one of following classes: \code{lba.ls}, \code{lba.ls.fe}, \code{lba.ls.logit}, \code{lba.mle}, \code{lba.mle.fe}, \code{lba.mle.logit}}
 \item{\dots}{ Further arguments (required by generic).}
}

\value{
  The \code{goodnessfit} function of the method \code{lba.mle}, \code{lba.mle.fe} and \code{lba.mle.logit} returns a list with the slots:
      
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

 The \code{goodnessfit} function of the method \code{lba.ls}, \code{lba.ls.fe} and \code{lba.ls.logit} returns a list with the slots:

 \item{dfdb     }{Degrees of freedom of the base model}
 \item{dfd      }{Degrees of freedom of the full model}
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
}

\references{
   Agresti, Alan. 2002. \emph{Categorical Data Analysis, second edition}. Hoboken: John Wiley and Sons.

   van der Ark, A. L. 1999. \emph{Contributions to Latent Budget Analysis, a tool for the analysis of compositional data}. Ph.D. Thesis University of Utrecht.
}

\note{For a detailed and complete discussion about goodness of fit results for latent budget analysis, see van der Ark 1999.
}

\seealso{
 \code{\link{summary.goodnessfit.lba.ls}},  \code{\link{summary.goodnessfit.lba.mle}},\code{\link{lba}}
}

\examples{
data('votB')

# Using LS method (default) without constraint
# K = 2
ex1 <- lba(parties ~ city,
           votB,
           K = 2)

gx1 <- goodnessfit(ex1)
gx1

# Using MLE method without constraint
# K = 2
exm <- lba(parties ~ city,
           votB,
           K = 2,
           method='mle')

gxm <- goodnessfit(exm)
gxm

# Using LS method (default) with LOGIT constrain
data('housing')

# Make cross-table to matrix design.
tbh <- xtabs(value ~ Influence + Housing, housing)

Xis <- model.matrix(~ Housing*Influence,
                    tbh,
                    contrasts=list(Housing='contr.sum',
                                   Influence='contr.sum'))

tby <- xtabs(value ~ Satisfaction + Contact, housing)

Yis <- model.matrix(~ Satisfaction*Contact,
                    tby,
                    contrasts=list(Satisfaction='contr.sum',
                                   Contact='contr.sum'))[,-1]

S <- 12
T <- 5

tabs <- xtabs(value ~ interaction(Housing,
                                  Influence) + interaction(Satisfaction,
                                  Contact),
              housing)
\dontrun{
ex2 <- lba(tabs,
           K = 2,
           logitA = Xis,
           logitB = Yis,
           S = S,
           T = T,
           trace.lba=FALSE) 

gex2 <- goodnessfit(ex2)
gex2
}
}
