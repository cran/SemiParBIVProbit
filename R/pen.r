pen <- function(params, qu.mag, sp, VC){


##############################################################
    ma1 <- matrix(0,VC$gp1,VC$gp1)   
    if(VC$l.sp1 == 0) EQ1P <- adiag(ma1)
    
    if(VC$l.sp1 != 0){
    ind <- 1:VC$l.sp1
    S1 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S1 <- do.call(adiag, lapply(S1, unlist))
    EQ1P <- adiag(ma1, S1)
                   } 
                   
##############################################################
    
    
    
##############################################################    
    ma2 <- matrix(0,VC$gp2,VC$gp2) 
    if(VC$l.sp2 == 0) EQ2P <- adiag(ma2)    
    
    if(VC$l.sp2 != 0){
    ind <- (VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp2)
    S2 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S2 <- do.call(adiag, lapply(S2, unlist))
    EQ2P <- adiag(ma2, S2)
                   }
                            
##############################################################
    
    
    
##############################################################    
    if(is.null(VC$gp3)) EQ3P <- 0
    
    if(!is.null(VC$gp3)){
    
    ma3 <- matrix(0,VC$gp3,VC$gp3) 
    
    if(VC$l.sp3 != 0){
    ind <- (VC$l.sp1 + VC$l.sp2 + 1):(VC$l.sp1 + VC$l.sp2 + VC$l.sp3)
    S3 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S3 <- do.call(adiag, lapply(S3, unlist))
    EQ3P <- adiag(ma3, S3)
                   }
                   
    if(VC$l.sp3 == 0) EQ3P <- adiag(ma3)      
    
    }
##############################################################
    

    S.h <- adiag(EQ1P, EQ2P, EQ3P)
    
    S.h1 <- 0.5*crossprod(params,S.h)%*%params
    S.h2 <- S.h%*%params
    
    list(S.h = S.h, S.h1 = S.h1, S.h2 = S.h2)
   
         }





















