goodnessfit.lba.ls <- function(object,...){

  stopifnot(inherits(object,'lba'))

  base <- update(object,
                 K = 1,
                 trace.lba=F)

  ifelse(inherits(object,'lba.formula'),
         N <- object$tab,
         N <-  eval(getCall(object)$obj))

  N[N==0] <- 1e-4
  I <- nrow(object[[4]])
  J <- nrow(object[[5]])
  K <- ncol(object[[4]])
  P <- object[[1]]
  pij <- object[[2]]

  row.weights <- eval(getCall(object)$row.weights)
  col.weights <- eval(getCall(object)$col.weights)

  if(is.null(row.weights)){
    vI <- rep(1,I)
    V  <- vI * diag(I)
  } else {
    vI <- row.weights
    V <- vI * diag(I)
  }

  if(is.null(col.weights)){
    wi <- rep(1,J)
    W  <- wi * diag(J)
  } else {
    wi <- col.weights
    W <- wi * diag(J)
  }

  pip <- rowSums(N)/sum(N)

  dfd <- (I-K)*(J-K) #degrees of freedom identified solutions 

  #K = 1  baseline model
  dfdb <- (I-1)*(J-1)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # GFS DISTRIBUTION FREE
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #RSS = residual sum of squares   
  #K=1 baseline    
  RSS1 <- sum((V%*%(base[[2]] - P)%*%W)^2)

  #K=K
  RSS <-  sum((V%*%(pij - P)%*%W)^2)

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

  class(res) <- 'goodnessfit.lba.ls'
  return(res)

} 
