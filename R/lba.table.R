lba.table <- function(obj,
                      A           = NULL, # mixing parameters
                      B           = NULL, # latent components
                      K           =  1L,  # integer 
                      cA          = NULL, # position of the constraints mixing parameters
                      cB          = NULL, # position of the constraints latent components 
                      logitA      = NULL, # design matrix for row covariates IxS 
                      logitB      = NULL, # design matrix for row covariates JxT 
                      omsk        = NULL, # matrix of logit parameters SxK
                      psitk       = NULL, # matrix of logit parameters TxK
                      S           = NULL, # integer
                      T           = NULL, # integer
                      row.weights = NULL, # row weights matrix
                      col.weights = NULL, # column weights matrix
                      tolG        = 1e-10, 
                      tolA        = 1e-05, 
                      tolB        = 1e-05,
                      itmax.unide = 1e3, # With and without constraint
                      itmax.ide   = 1e3,
                      trace.lba   = TRUE,# only when K > 3
                      toltype     = 'all',
                      method      = c("ls", 
                                      "mle"),
                      what        = c('inner',
                                      'outer'),
                      ...) 
{

  if(length(dim(obj)) == 3L){

    stop('Your table should be one-dimensional!')

  }

  switch(match.arg(what),
         inner = what <- 'inner',
         outer = what <- 'outer')

  switch(match.arg(method),
         ls = method <- 'ls',
         mle = method <- 'mle')

  if(is.null(cA) & is.null(cB) & is.null(logitA) & is.null(logitB)){

    class(obj)  <- method

    result <- lba(obj,
                  A           =  A,           
                  B           =  B,
                  K           =  K,
                  row.weights =  row.weights, 
                  col.weights =  col.weights, 
                  tolG        =  tolG,        
                  tolA        =  tolA,        
                  tolB        =  tolB,        
                  itmax.unide =  itmax.unide,
                  itmax.ide   =  itmax.ide,
                  trace.lba   =  trace.lba,  
                  toltype     =  toltype,
                  what        =  what,
                  ...) 

  } else 

    if((!is.null(cA) | !is.null(cB)) & is.null(logitA) & is.null(logitB)){

      class(obj) <- paste(method,
                          'fe',
                          sep='.')

      result <- lba(obj,
                    A           =  A,           
                    B           =  B,
                    K           =  K,
                    cA          =  cA,          
                    cB          =  cB,           
                    row.weights =  row.weights, 
                    col.weights =  col.weights, 
                    tolG        =  tolG,        
                    tolA        =  tolA,        
                    tolB        =  tolB,        
                    itmax.ide   =  itmax.ide,
                    trace.lba   =  trace.lba,  
                    toltype     =  toltype,
                    what        =  what,
                    ...)

    } else {

      class(obj) <- paste(method,
                          'logit',
                          sep='.')

      result <- lba(obj,
                    A           =  A,           
                    B           =  B,
                    K           =  K,
                    cA          =  cA,
                    cB          =  cB,
                    logitA      =  logitA,      
                    logitB      =  logitB,      
                    omsk        =  omsk,        
                    psitk       =  psitk,       
                    S           =  S,           
                    T           =  T,           
                    row.weights =  row.weights, 
                    col.weights =  col.weights, 
                    tolG        =  tolG,        
                    tolA        =  tolA,        
                    tolB        =  tolB,        
                    itmax.ide   =  itmax.ide,
                    trace.lba   =  trace.lba,  
                    toltype     =  toltype,
                    what        =  what,
                    ...)

    }

  n_dim <- length(result$pk)-1

  if(n_dim == 1){

    class(result) <- c('lba.1d',
                       class(result),
                       'lba.table',
                       'lba')
  }

  if(n_dim == 2){

    class(result) <- c('lba.2d',
                       class(result),
                       'lba.table',
                       'lba')
  }

  if(n_dim >= 3){

    class(result) <- c('lba.3d',
                       class(result),
                       'lba.table',
                       'lba')
  } else {

    class(result) <- c(class(result),
                       'lba.table',
                       'lba') 

  }

  cl <- match.call()

  result$call <- cl
  result$what <- what

  result
}
