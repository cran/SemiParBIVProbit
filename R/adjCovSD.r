adjCovSD <- function(object,design){

# this does not work with power links and two-parameter copulas


    Ainv <- object$Vb 
    estfun <- cbind( c(object$fit$dl.dbe1)*object$X1, c(object$fit$dl.dbe2)*object$X2, object$fit$dl.drho )
    
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
  
  