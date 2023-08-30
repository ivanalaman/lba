goodnessfit.lba.mle.fe <- function(object,...){

  stopifnot(inherits(object,'lba'))

  base <- update(object,
                 cA=NULL,
                 cB=NULL,
                 trace.lba=F,
                 K = 1)

  ifelse(inherits(object,'lba.formula'),
         N <- object$tab,
         N <-  eval(getCall(object)$obj))

  N[N==0] <- 1e-4
  G2b <- 2 * sum(N * log(N/(base[[2]] * rowSums(N))))

  Nexpected <-  base[[2]] * rowSums(N)
  chi2b <- sum(((N - Nexpected )^2)/Nexpected) 

  I <- nrow(object[[4]])
  J <- nrow(object[[5]])
  K <- ncol(object[[4]])
  P <- object[[1]]
  pij <- object[[2]]

  pip <- rowSums(N)/sum(N)

  G2 <- 2 * sum(N * log(N/(pij * rowSums(N))))

  pexpected <-  pij * rowSums(N)
  chi2 <- sum(((N - pexpected)^2)/pexpected) 

  #++++ calculus of dfd
  recu_cA = eval(getCall(object)$cA) # recuperando a matrix de restricao
  recu_cB = eval(getCall(object)$cB) # recuperando a matrix de restricao

  if(!is.null(recu_cA)){
       au1 = as.vector(recu_cA) 

       #+++++ONLY EQUALITY++++###
       au1e = au1[!is.na(au1) & au1 > 1]
       au2e = table(au1e)

       ifelse(length(au2e) > 0,
              au3e <- sum(sapply(au2e,function(x) x-1)),
              au3e <- 0)

       #+++++ONLY FIXED+++++###
       au1f = au1[!is.na(au1) & au1 <= 1]
       au2f = length(au1f)

       dfA <- au3e + au2f
         } else {
           dfA <- 0
         }

  if(!is.null(recu_cB)){
    bu1 = as.vector(recu_cB) 

    bu1e = bu1[!is.na(bu1) & bu1 > 1]
    bu2e = table(bu1e)

    ifelse(length(bu2e) > 0,
           bu3e <- sum(sapply(bu2e,function(x) x-1)),
           bu3e <- 0)

    bu1f = bu1[!is.na(bu1) & bu1 <= 1]
    bu2f = length(bu1f)

    dfB <- bu3e + bu2f
  } else {
    dfB <- 0
  }

  dfd <- I*(J-1)-I*(K-1)-K*(J-1)+dfA+dfB #degrees of freedom identified solutions 

  #K = 1  baseline model
  dfdb <- (I-1)*(J-1)

  #K = K  p-value
  proG <- 1-pchisq(G2, 
                   df = dfd) #p-value of G2

  prochi <- 1-pchisq(chi2, 
                     df = dfd) #p-value of chi2  

  #K = 1
  proG1 <- 1-pchisq(G2b, 
                    df = dfdb) #p-value of G2

  prochi1 <- 1-pchisq(chi2b, 
                      df = dfdb) #p-value of chi2 

  #improvement
  impG_mle <- abs(G2b - G2) 

  # improvement per budget page 164
  impPB_mle <- G2b/min(I,J)

  #average improvement per degree of freedom
  impDF_mle <- G2b/((I-1)*(J-1))                     

  #==============================================================================
  #          Low value  GFS
  #==============================================================================
  AICC <- G2 - 2*dfd

  BICC <- G2 - dfd*log(I*J)

  CAIC <- G2 - dfd*(log(I*J) + 1)

  AICb <- G2b - 2*dfdb

  BICb <- G2b - dfdb*log(I*J)

  CAICb <- G2b - dfdb*(log(I*J) + 1)
  #==============================================================================

  #==============================================================================
  #         Incremental GFS
  #==============================================================================
  delta1 <- (G2b - G2)/G2b
  delta2 <- (G2b - G2)/(G2b-dfd)
  rho1 <- (G2b/dfdb - G2/dfd)/(G2b/dfdb)                 
  rho2 <- (G2b/dfdb - G2/dfd)/(G2b/dfdb - 1)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # GFS DISTRIBUTION FREE
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #RSS = residual sum of squares   
  #K=1 baseline
  RSS1 <- sum((base[[2]] - P)^2)

  #K=K
  RSS <- sum((pij - P)^2)

  #improvement
  impRSS <- RSS1 - RSS

  # improvement per budget page 164
  impPB <- RSS1/min(I,J)

  #average improvement per degree of freedom
  impDF <- RSS1/((I-1)*(J-1))

  #Index of Dissimilarity; D
  #K = 1 baseline
  D1 <- 0.5*sum(pip*rowSums(abs(base[[2]] - P)))

  #K = K
  D <- 0.5*sum(pip*rowSums(abs(pij - P)))

  #Proportion of correctly classified data; 1 - D
  #K = 1
  pccb <- 1-D1

  #K = K
  pcc <- 1 - D 

  #improvement
  impD <- (D1 - D)/(pccb)

  #improvement of Proportion of correctly classified data per budget
  impPCCB <- (pccb)/min(I,J)

  #average improvement of Proportion of correctly classified data per degree of freedom
  AimpPCCDF <- (pccb)/((I-1)*(J-1))  

  # m.a.d; mean angular deviation

  #K = 1; baseline
  mad1 <- sum(acos(rowSums(base[[2]]*P)/sqrt(rowSums(P^2)*rowSums(base[[2]]^2))))/I

  #K = K
  madk <- sum(acos(rowSums(pij*P)/sqrt(rowSums(P^2)*rowSums(pij^2))))/I

  #improvement
  impMad <- (mad1 - madk)/mad1

  # improvement per budget
  impPBsat <- mad1/min(I,J)   
  #in relation to SAT in which case is not zero

  #average improvement per degree of freedom
  impDFsat <- mad1/((I-1)*(J-1))

  res <- list(dfdb,
              dfd,
              G2b,
              G2,
              chi2b,
              chi2,
              proG1,
              proG,
              prochi1,
              prochi,
              AICb,
              AICC,
              BICb,
              BICC,
              CAICb,
              CAIC,
              delta1,
              delta2,
              rho1,
              rho2,
              RSS1,
              RSS,
              impRSS,
              impPB,
              impDF,
              D1,
              D,
              pccb,
              pcc,
              impD,
              impPCCB,
              AimpPCCDF,
              mad1,
              madk,
              impMad,
              impPBsat,
              impDFsat)

  names(res) <- c('dfdb',
                  'dfd',
                  'G2b',
                  'G2',
                  'chi2b',
                  'chi2',
                  'proG1',
                  'proG',
                  'prochi1',
                  'prochi',
                  'AICb',
                  'AICC',
                  'BICb',
                  'BICC',
                  'CAICb',
                  'CAIC',
                  'delta1',
                  'delta2',
                  'rho1',
                  'rho2', 
                  'RSS1',
                  'RSS',
                  'impRSS',
                  'impPB',
                  'impDF',
                  'D1', 
                  'D',
                  'pccb',
                  'pcc',
                  'impD',
                  'impPCCB',
                  'AimpPCCDF',
                  'mad1',
                  'madk',
                  'impMad',
                  'impPBsat',
                  'impDFsat') 

  class(res) <- c('goodnessfit.lba.mle.fe',
                  'goodnessfit',
                  'list')
  return(res) 
}
