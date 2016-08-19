bprobgHsContUniv <- function(params, respvec, VC, ps, AT = FALSE){

epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)

weights <- VC$weights

#########################################################################

if(VC$Cont == "NO"){

X2 <- VC$X2
X3 <- VC$X3
y2 <- respvec$y2

if(VC$ccss == "yes"){weights <- weights[VC$inde]} #; y2 <- y2[VC$inde]; X2 <- X2[VC$inde,]; X3 <- X3[VC$inde,] } 


eta2 <- X2%*%params[1:VC$X2.d2]
eta2 <- eta.tr(eta2, VC$margins[2])



if( !(VC$margins[2] %in% c(VC$m1d)) ){ 


if(is.null(VC$X3))  sigma2.st <- params[(VC$X2.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- X3%*%params[(VC$X2.d2+1):(VC$X2.d2+VC$X3.d2)] 
              
}              
              
if(VC$margins[2] %in% VC$m1d) sigma2.st <- 0 # as there is no sigma              
              
                  
    sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 


if(VC$margins[2] %in% VC$m2)            dHs  <-      distrHs(y2, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = TRUE)
if(VC$margins[2] %in% c(VC$m1d,VC$m2d)) dHs  <- distrHsDiscr(y2, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = TRUE, y2m = VC$y2m)


}

#########################################################################
#########################################################################


if(VC$Cont == "YES" && respvec$univ == 2){ # interested in first equation


eta2 <- VC$X1%*%params[1:VC$X1.d2] # this is eta1 but changed to eta2 for convenience
eta2 <- eta.tr(eta2, VC$margins[1])


if( !(VC$margins[1] %in% c(VC$m1d)) ){ 

if(is.null(VC$X3))  sigma2.st <- params[(VC$X1.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- VC$X3%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]

}

if(VC$margins[1] %in% VC$m1d) sigma2.st <- 0



    sstr1 <- esp.tr(sigma2.st, VC$margins[1])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 


if(VC$margins[1] %in% VC$m2)            dHs  <-      distrHs(respvec$y1, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = TRUE)
if(VC$margins[1] %in% c(VC$m1d,VC$m2d)) dHs  <- distrHsDiscr(respvec$y1, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = TRUE, y2m = VC$y1m)

}


if(VC$Cont == "YES" && respvec$univ == 3){ # interested in second equation


eta2 <- VC$X2%*%params[1:VC$X2.d2]
eta2 <- eta.tr(eta2, VC$margins[2])


if( !(VC$margins[2] %in% c(VC$m1d)) ){ 

if(is.null(VC$X3))  sigma2.st <- params[(VC$X2.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- VC$X4%*%params[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]

}


if(VC$margins[2] %in% VC$m1d) sigma2.st <- 0


    sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 

if(VC$margins[2] %in% VC$m2)            dHs <-      distrHs(respvec$y2, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = TRUE)
if(VC$margins[2] %in% c(VC$m1d,VC$m2d)) dHs <- distrHsDiscr(respvec$y2, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = TRUE, y2m = VC$y2m)


}

#########################################################################
#########################################################################


pdf2                         <- dHs$pdf2
derpdf2.dereta2              <- dHs$derpdf2.dereta2 
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st  
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2         
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  
  
    
########################################################################################################

l.par <- weights*log(pdf2)
res   <- -sum(l.par)
  
########################################################################################################
 
dl.dbe       <- -weights*( derpdf2.dereta2/pdf2 )
dl.dsigma.st <- -weights*( derpdf2.dersigma2.st/pdf2 )
                     
d2l.be.be        <- -weights*( (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 )
d2l.sigma.sigma  <- -weights*( (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 )
d2l.be.sigma     <- -weights*( (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 )
     
########################################################################################################
     

if(VC$Cont == "NO"){


if( !(VC$margins[2] %in% c(VC$m1d)) ){ ##

if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*X2 ) ,
            sum( dl.dsigma.st ) )
                
  be.be    <- crossprod(X2*c(d2l.be.be),X2)
  be.sigma <- t(t(rowSums(t(X2*c(d2l.be.sigma))))) 

  H <- rbind( cbind( be.be      , be.sigma             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma) )  ) 
                            
                    }


if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*X2        ) ,
            colSums( c(dl.dsigma.st)*X3  ) )
                
  be.be    <- crossprod(X2*c(d2l.be.be),X2)
  be.sigma <- crossprod(X2*c(d2l.be.sigma),X3)
  si.sigma <- crossprod(X3*c(d2l.sigma.sigma),X3)  

  H <- rbind( cbind( be.be      , be.sigma  ), 
              cbind( t(be.sigma), si.sigma  )  ) 
              
                      }

} ##


if( VC$margins[2] %in% c(VC$m1d) ){

  G <- c( colSums( c(dl.dbe)*X2 ) )          
  H <- crossprod(X2*c(d2l.be.be),X2)

}


}



if(VC$Cont == "YES" && respvec$univ == 2){



if( !(VC$margins[1] %in% c(VC$m1d)) ){ ##

if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X1 ) ,
            sum( dl.dsigma.st )
          )
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be.sigma))))) 

  H <- rbind( cbind( be.be      , be.sigma             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma) )  
            ) 
                    }


if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X1        ) ,
            colSums( c(dl.dsigma.st)*VC$X3  )
          )
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- crossprod(VC$X1*c(d2l.be.sigma),VC$X3)
  si.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)  

  H <- rbind( cbind( be.be      , be.sigma  ), 
              cbind( t(be.sigma), si.sigma  )  
            )  
                     }
                                         
                     

} ##


if( VC$margins[1] %in% c(VC$m1d) ){

  G <- c( colSums( c(dl.dbe)*VC$X1 ) )          
  H <- crossprod(VC$X1*c(d2l.be.be),VC$X1)

}




}





if(VC$Cont == "YES" && respvec$univ == 3){


if( !(VC$margins[2] %in% c(VC$m1d)) ){ ##

if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2 ) ,
            sum( dl.dsigma.st )
          )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be.sigma))))) 

  H <- rbind( cbind( be.be      , be.sigma             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma) )  
            ) 
                     }


if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2        ) ,
            colSums( c(dl.dsigma.st)*VC$X4  )
          )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- crossprod(VC$X2*c(d2l.be.sigma),VC$X4)
  si.sigma <- crossprod(VC$X4*c(d2l.sigma.sigma),VC$X4)  

  H <- rbind( cbind( be.be      , be.sigma  ), 
              cbind( t(be.sigma), si.sigma  )  
            )  
                      }




}



if( VC$margins[2] %in% c(VC$m1d) ){ # this useful for later

  G <- c( colSums( c(dl.dbe)*VC$X2 )  )         
  H <- crossprod(VC$X2*c(d2l.be.be),VC$X2)

}




}




if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
  S.h <- ps$S.h  

  if( length(S.h) != 1 ){
  
  S.h1 <- 0.5*crossprod(params,S.h)%*%params
  S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   

  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h 
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2) 
  

list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.par=l.par, ps = ps, sigma2.st = sigma2.st,
     BivD=VC$BivD, eta2 = eta2, sigma2 = sigma2, nu = NULL)      


}


