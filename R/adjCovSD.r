adjCovSD <- function(object, design){

if(object$PL!="P") stop("This correction does not currently work for models with asymmetric links.")

    Ainv <- object$Vb 
    
    if(is.null(object$X3) )  mul <- 1
    if(!is.null(object$X3) ) mul <- object$X3[object$good,]
    
    estfun <- cbind( c(object$fit$dl.dbe1)*object$X1[object$good,], c(object$fit$dl.dbe2)*object$X2[object$good,], object$fit$dl.drho*mul )
    
    if (inherits(design,"survey.design2"))
      covsan <- svyrecvar(estfun%*%Ainv,design$cluster,design$strata,design$fpc,postStrata=design$postStrata)
    else if (inherits(design, "twophase"))
      covsan <- twophasevar(estfun%*%Ainv, design)
    else if (inherits(design, "twophase2"))
      covsan <- twophase2var(estfun%*%Ainv, design)
    else
      covsan <- svyCprod(estfun%*%Ainv,design$strata,design$cluster[[1]],design$fpc, design$nPSU,
                  design$certainty,design$postStrata)
                 
  object$Vb <- covsan
  object              
                                             
  }
  
  