pen <- function(params, qu.mag, sp, VC, eq1 = "yes"){

cont2par <- c("N","GU","rGU","LO","LN","WEI","WEI2","iG","GA") 
cont3par <- c("DAGUM") 

l.sp1 <- VC$l.sp1 
l.sp2 <- VC$l.sp2 
l.sp3 <- VC$l.sp3 
l.sp4 <- VC$l.sp4 
l.sp5 <- VC$l.sp5 
l.sp6 <- VC$l.sp6 


  if(eq1 == "no"){   
  
    if(VC$margins[2] %in% cont2par)  l.sp1 <- l.sp4 <- l.sp5 <- l.sp6 <- 0 
    if(VC$margins[2] %in% cont3par)  l.sp1 <- l.sp5 <- l.sp6 <- 0 
 
  } 





##############################################################
    ma1 <- matrix(0,VC$gp1,VC$gp1)   
    if(l.sp1 == 0) EQ1P <- adiag(ma1)
    
    if(l.sp1 != 0){
    ind <- 1:l.sp1
    S1 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S1 <- do.call(adiag, lapply(S1, unlist))
    EQ1P <- adiag(ma1, S1)
                   } 
                   
##############################################################
    
##############################################################    
    ma2 <- matrix(0,VC$gp2,VC$gp2) 
    if(l.sp2 == 0) EQ2P <- adiag(ma2)    
    
    if(l.sp2 != 0){
    ind <- (l.sp1 + 1):(l.sp1 + l.sp2)
    S2 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S2 <- do.call(adiag, lapply(S2, unlist))
    EQ2P <- adiag(ma2, S2)
                   }
                            
##############################################################
    
##############################################################
    
if(!is.null(VC$gp3)){

 EQ4P <- EQ5P <- EQ6P <- NULL 

    ma3 <- matrix(0,VC$gp3,VC$gp3) 
    if(l.sp3 == 0) EQ3P <- adiag(ma3)    
    
    if(l.sp3 != 0){
    ind <- (l.sp1 + l.sp2 + 1):(l.sp1 + l.sp2 + l.sp3)
    S3 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S3 <- do.call(adiag, lapply(S3, unlist))
    EQ3P <- adiag(ma3, S3)
                   }
        
    if(!is.null(VC$gp4)){
    
    ma4 <- matrix(0,VC$gp4,VC$gp4) 
    if(l.sp4 == 0) EQ4P <- adiag(ma4)    
    
    if(l.sp4 != 0){
    ind <- (l.sp1 + l.sp2 + l.sp3 + 1):(l.sp1 + l.sp2 + l.sp3 + l.sp4)
    S4 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S4 <- do.call(adiag, lapply(S4, unlist))
    EQ4P <- adiag(ma4, S4)
                   }
                        }
                        
	    if(!is.null(VC$gp5)){
	    
	    ma5 <- matrix(0,VC$gp5,VC$gp5) 
	    if(l.sp5 == 0) EQ5P <- adiag(ma5)    
	    
	    if(l.sp5 != 0){
	    ind <- (l.sp1 + l.sp2 + l.sp3 + l.sp4 + 1):(l.sp1 + l.sp2 + l.sp3 + l.sp4 + l.sp5)
	    S5 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    S5 <- do.call(adiag, lapply(S5, unlist))
	    EQ5P <- adiag(ma5, S5)
	                     } 
               
                                }
            if(!is.null(VC$gp6)){
	    
	    ma6 <- matrix(0,VC$gp6,VC$gp6) 
	    if(l.sp6 == 0) EQ6P <- adiag(ma6)    
	    
	    if(l.sp6 != 0){
	    ind <- (l.sp1 + l.sp2 + l.sp3 + l.sp4 + l.sp5 + 1):(l.sp1 + l.sp2 + l.sp3 + l.sp4 + l.sp5 + l.sp6)
	    S6 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    S6 <- do.call(adiag, lapply(S6, unlist))
	    EQ6P <- adiag(ma6, S6)
	                     } 
               
                                }  
                                
                                
}else{

    if(VC$Model == "BPO0")                               {EQ3P <- EQ4P <- EQ5P <- EQ6P <- NULL                 }
    if(VC$margins[2] == "probit" && VC$Model != "BPO0")  {EQ3P <- 0; EQ4P <- NULL; EQ5P <- NULL; EQ6P <- NULL  }
    if(VC$margins[2] %in% cont2par)                      {EQ3P <- 0; EQ4P <- 0;    EQ5P <- NULL; EQ6P <- NULL  }
    if(VC$margins[2] %in% cont3par)                      {EQ3P <- 0; EQ4P <- 0;    EQ5P <- 0;    EQ6P <- NULL  }



}  
    
    
    if( eq1 == "no" && VC$margins[2] %in% cont2par) {EQ4P <- NULL; S.h <- adiag(EQ2P, EQ3P, EQ4P, EQ5P, EQ6P)}  
    if( eq1 == "no" && VC$margins[2] %in% cont3par) {EQ5P <- NULL; S.h <- adiag(EQ2P, EQ3P, EQ4P, EQ5P, EQ6P)}  
    
    if( eq1 == "yes") S.h <- adiag(EQ1P, EQ2P, EQ3P, EQ4P, EQ5P, EQ6P)
    
    
    S.h1 <- 0.5*crossprod(params,S.h)%*%params
    S.h2 <- S.h%*%params
    
    list(S.h = S.h, S.h1 = S.h1, S.h2 = S.h2)
   
         }





















