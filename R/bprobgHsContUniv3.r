bprobgHsContUniv3 <- function(params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){


eta2 <- VC$X2%*%params[1:VC$X2.d2]

epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)

if(is.null(VC$X3))  sigma2.st <- params[(VC$X2.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- VC$X3%*%params[(VC$X2.d2+1):(VC$X2.d2+VC$X3.d2)]

if(is.null(VC$X4))  nu.st <- params[(VC$X2.d2 + 2)] 
if(!is.null(VC$X4)) nu.st <- VC$X4%*%params[(VC$X2.d2+VC$X3.d2+1):(VC$X2.d2+VC$X3.d2+VC$X4.d2)]

#sigma2.st <- ifelse(sigma2.st > 4, 4, sigma2.st) # new 

sigma2 <- exp(sigma2.st) + epsilon


  if(VC$margins[2] == "DAGUM") nu <- exp(nu.st) + epsilon



  eta2 <- ifelse( eta2 > 600, 600, eta2 )
  eta2 <- ifelse( eta2 < -17, -17, eta2 )  
  
dHs  <- distrHs(respvec$y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2=VC$margins[2], naive = TRUE)
  
pdf2                         <- dHs$pdf2
derpdf2.dereta2              <- dHs$derpdf2.dereta2 
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st  
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2         
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  


der2pdf2.dereta2dernu.st    <- dHs$der2pdf2.dereta2dernu.st   
der2pdf2.dersigma2.stdernu.st <- dHs$der2pdf2.sigma2.st2dernu.st
derpdf2.dernu.st            <- dHs$derpdf2.dernu.st           
der2pdf2.dernu.st2          <- dHs$der2pdf2.dernu.st2         
    
########################################################################################################

l.par <- VC$weights*log(pdf2)
res   <- -sum(l.par)
  
########################################################################################################
 
dl.dbe       <- -VC$weights*( derpdf2.dereta2/pdf2 )
dl.dsigma.st <- -VC$weights*( derpdf2.dersigma2.st/pdf2 )
dl.dnu.st    <- -VC$weights*( derpdf2.dernu.st/pdf2 )   
          
          
d2l.be.be        <- -VC$weights*( (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 )
d2l.sigma.sigma  <- -VC$weights*( (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 )
d2l.be.sigma     <- -VC$weights*( (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 )
d2l.be.nu        <- -VC$weights*( (der2pdf2.dereta2dernu.st*pdf2 - derpdf2.dereta2*derpdf2.dernu.st)/pdf2^2 )
     
d2l.sigma.nu     <- -VC$weights*( (der2pdf2.dersigma2.stdernu.st*pdf2-(derpdf2.dersigma2.st*derpdf2.dernu.st))/pdf2^2 )
  d2l.nu.nu      <- -VC$weights*( (der2pdf2.dernu.st2*pdf2-(derpdf2.dernu.st)^2)/pdf2^2 )
     
     
     
########################################################################################################
     

if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2 ) ,
            sum( dl.dsigma.st ), sum( dl.dnu.st ) 
          )
                

  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be.sigma))))) 
  be.nu    <- t(t(rowSums(t(VC$X2*c(d2l.be.nu))))) 


  H <- rbind( cbind( be.be      , be.sigma,             be.nu            ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma), sum(d2l.sigma.nu)                  )  ,
              cbind( t(be.nu),    sum(d2l.sigma.nu),    sum(d2l.nu.nu)              )
              
            ) 
  
}




if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2        ) ,
            colSums( c(dl.dsigma.st)*VC$X3  ),
            colSums( c(dl.dnu.st)*VC$X4  )
          )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- crossprod(VC$X2*c(d2l.be.sigma),VC$X3)
  be.nu    <- crossprod(VC$X2*c(d2l.be.nu),VC$X4)
  
  si.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)  
  si.nu    <- crossprod(VC$X4*c(d2l.nu.nu),VC$X4)    
  
  sa.nu    <- crossprod(VC$X3*c(d2l.sigma.nu),VC$X4)    
  
  
  

  H <- rbind( cbind( be.be      , be.sigma , be.nu ), 
              cbind( t(be.sigma), si.sigma,  sa.nu ),
              cbind( t(be.nu),    t(sa.nu),  si.nu ) 

            )  
   
    
}



if( (VC$l.sp2==0 && VC$l.sp3==0 && VC$l.sp4==0) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, eq1 = "no")


if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2) 
  

         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps, sigma2.st = sigma2.st, nu.st = nu.st,
              eta2=eta2, 
dl.dbe          = dl.dbe,            
dl.dsigma.st    = dl.dsigma.st,
d2l.be.be       = d2l.be.be,
d2l.sigma.sigma = d2l.sigma.sigma,
d2l.be.sigma    = d2l.be.sigma,
              BivD=VC$BivD)      

}




     























