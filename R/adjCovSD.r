adjCovSD <- function(x, design){



cont2par <- c("N","GU","rGU","LO","LN","WEI","iG","GA","iGA") 

Ainv <- x$Vb 
    
if( x$margins[2]=="probit" && x$Model != "BPO0"){

if(is.null(x$X3) )  mul <- 1
if(!is.null(x$X3) ) mul <- x$X3[x$good,]

estfun <- cbind( c(x$fit$dl.dbe1)*x$X1[x$good,], c(x$fit$dl.dbe2)*x$X2[x$good,], c(x$fit$dl.drho)*mul )

}

if( x$Model == "BPO0" ){

estfun <- cbind( c(x$fit$dl.dbe1)*x$X1[x$good,], c(x$fit$dl.dbe2)*x$X2[x$good,] )

}



if( x$margins[2] %in% cont2par ){

if( !is.null(x$X3) && !is.null(x$X4) ) mul1 <- x$X3; mul2 <- x$X4 
if(  is.null(x$X3) &&  is.null(x$X4) ) mul1 <- mul2 <- 1 
                                       
estfun <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma.st)*mul1,
                 c(x$fit$dl.dteta.st)*mul2       )                                           
}

    
    

    if (inherits(design,"survey.design2"))
      covsan <- svyrecvar(estfun%*%Ainv,design$cluster,design$strata,design$fpc,postStrata=design$postStrata)
    else if (inherits(design, "twophase"))
      covsan <- twophasevar(estfun%*%Ainv, design)
    else if (inherits(design, "twophase2"))
      covsan <- twophase2var(estfun%*%Ainv, design)
    else
      covsan <- svyCprod(estfun%*%Ainv,design$strata,design$cluster[[1]],design$fpc, design$nPSU,
                  design$certainty,design$postStrata)
                 
  x$Vb <- covsan
  
  rm(Ainv, estfun, covsan)
  
  x              
                                             
  }
  
  