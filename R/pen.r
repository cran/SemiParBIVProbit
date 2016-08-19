pen <- function(qu.mag, sp, VC, univ, l.splist){


##############################################################
    ma1 <- matrix(0,VC$gp1,VC$gp1)   
    if(l.splist$l.sp1 == 0) EQ1P <- adiag(ma1)
    
    if(l.splist$l.sp1 != 0){
    ind <- 1:l.splist$l.sp1
    S1 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S1 <- do.call(adiag, lapply(S1, unlist))
    EQ1P <- adiag(ma1, S1)
                   } 
                   
##############################################################
    
##############################################################    
    ma2 <- matrix(0,VC$gp2,VC$gp2) 
    if(l.splist$l.sp2 == 0) EQ2P <- adiag(ma2)    
    
    if(l.splist$l.sp2 != 0){
    ind <- (l.splist$l.sp1 + 1):(l.splist$l.sp1 + l.splist$l.sp2)
    S2 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S2 <- do.call(adiag, lapply(S2, unlist))
    EQ2P <- adiag(ma2, S2)
                   }
                            
##############################################################
    
##############################################################
    
if(!is.null(VC$gp3)){ # this starts after first two equations

 EQ4P <- EQ5P <- EQ6P <- EQ7P <- NULL 

    ma3 <- matrix(0,VC$gp3,VC$gp3) 
    if(l.splist$l.sp3 == 0) EQ3P <- adiag(ma3)    
    
    if(l.splist$l.sp3 != 0){
    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3)
    S3 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S3 <- do.call(adiag, lapply(S3, unlist))
    EQ3P <- adiag(ma3, S3)
                   }
        
    if(!is.null(VC$gp4)){
    
    ma4 <- matrix(0,VC$gp4,VC$gp4) 
    if(l.splist$l.sp4 == 0) EQ4P <- adiag(ma4)    
    
    if(l.splist$l.sp4 != 0){
    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4)
    S4 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S4 <- do.call(adiag, lapply(S4, unlist))
    EQ4P <- adiag(ma4, S4)
                   }
                        }
                        
	    if(!is.null(VC$gp5)){
	    
	    ma5 <- matrix(0,VC$gp5,VC$gp5) 
	    if(l.splist$l.sp5 == 0) EQ5P <- adiag(ma5)    
	    
	    if(l.splist$l.sp5 != 0){
	    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5)
	    S5 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    S5 <- do.call(adiag, lapply(S5, unlist))
	    EQ5P <- adiag(ma5, S5)
	                     } 
                                }
                                
            if(!is.null(VC$gp6)){
	    
	    ma6 <- matrix(0,VC$gp6,VC$gp6) 
	    if(l.splist$l.sp6 == 0) EQ6P <- adiag(ma6)    
	    
	    if(l.splist$l.sp6 != 0){
	    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6)
	    S6 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    S6 <- do.call(adiag, lapply(S6, unlist))
	    EQ6P <- adiag(ma6, S6)
	                     } 
                                }  
                                
            if(!is.null(VC$gp7)){
	    
	    ma7 <- matrix(0,VC$gp7,VC$gp7) 
	    if(l.splist$l.sp7 == 0) EQ7P <- adiag(ma7)    
	    
	    if(l.splist$l.sp7 != 0){
	    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + l.splist$l.sp7)
	    S7 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    S7 <- do.call(adiag, lapply(S7, unlist))
	    EQ7P <- adiag(ma7, S7)
	                     }               
                                }                                  
                                
                                
}else{

    if(VC$Model == "BPO0")                                                       {EQ3P <- EQ4P <- EQ5P <- EQ6P <- EQ7P <- NULL                 }
    if(VC$margins[1] %in% c(VC$bl)        && VC$Model != "BPO0")                 {EQ3P <- 0; EQ4P <- NULL; EQ5P <- NULL; EQ6P <- EQ7P <- NULL  }
    
    if(VC$margins[1] %in% c(VC$bl,VC$m1d) && VC$margins[2] %in% c(VC$m1d))       {EQ3P <- 0; EQ4P <- NULL; EQ5P <- NULL; EQ6P <- EQ7P <- NULL  }
    
    if(VC$margins[1] %in% c(VC$bl,VC$m1d) && VC$margins[2] %in% c(VC$m2,VC$m2d)) {EQ3P <- 0; EQ4P <- 0;    EQ5P <- NULL; EQ6P <- EQ7P <- NULL  }
    if(VC$margins[1] %in% c(VC$bl,VC$m1d) && VC$margins[2] %in% VC$m3)           {EQ3P <- 0; EQ4P <- 0;    EQ5P <- 0;    EQ6P <- EQ7P <- NULL  }
    

    if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% c(VC$m2,VC$m2d)) {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- NULL; EQ7P <- NULL  } 
    if(VC$margins[1] %in% VC$m3           && VC$margins[2] %in% VC$m3)           {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0;    EQ7P <- 0  } 
    if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% VC$m3)           {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0;    EQ7P <- NULL  }  
    if(VC$margins[1] %in% VC$m3           && VC$margins[2] %in% c(VC$m2,VC$m2d)) {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0;    EQ7P <- NULL  }  

    if(VC$Model == "B" && !is.null(VC$theta.fx))                                     {EQ3P <- EQ4P <- EQ5P <- EQ6P <- EQ7P <- NULL                 }
    # not efficient, it could be done above but done here for the moment


}  
    
    
    
    
    
    
    
    
    
    
    
if(VC$triv == FALSE){    ### TRIV
    
    
    
    if(univ == 0) S.h <- adiag(EQ1P, EQ2P, EQ3P, EQ4P, EQ5P, EQ6P, EQ7P)
    
    
    
    if(univ == 1){
    
       if(VC$margins[2] %in% c(VC$m1d))       {S.h <- adiag(EQ2P)}  
       if(VC$margins[2] %in% c(VC$m2,VC$m2d)) {S.h <- adiag(EQ2P, EQ3P)}  
       if(VC$margins[2] %in% VC$m3)           {S.h <- adiag(EQ2P, EQ3P, EQ4P)}  
    
                 }
    
    
    if(univ == 2){
    
    
       if(VC$margins[1] %in% c(VC$m1d) && VC$margins[2] %in% c(VC$m2,VC$m3,VC$m2d,VC$m1d)) {S.h <- adiag(EQ1P)        } 
       
       if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% c(VC$m2,VC$m2d)) {S.h <- adiag(EQ1P, EQ3P)        } 
       if(VC$margins[1] %in% VC$m3           && VC$margins[2] %in% c(VC$m2,VC$m2d)) {S.h <- adiag(EQ1P, EQ3P, EQ5P)  } 
       if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% VC$m3) {S.h <- adiag(EQ1P, EQ3P)        }  
       if(VC$margins[1] %in% VC$m3           && VC$margins[2] %in% VC$m3) {S.h <- adiag(EQ1P, EQ3P, EQ5P)  }     
    
    }
    
    
    if(univ == 3){ 
    
       if(VC$margins[1] %in% c(VC$m1d)       && VC$margins[2] %in% c(VC$m1d))       {S.h <- adiag(EQ2P)              }     
       if(VC$margins[1] %in% c(VC$m1d)       && VC$margins[2] %in% c(VC$m2,VC$m2d)) {S.h <- adiag(EQ2P, EQ3P)        } 
       if(VC$margins[1] %in% c(VC$m1d)       && VC$margins[2] %in% VC$m3)           {S.h <- adiag(EQ2P, EQ3P, EQ4P)  }  
       
    
       if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% c(VC$m2,VC$m2d)) {S.h <- adiag(EQ2P, EQ4P)        } 
       if(VC$margins[1] %in% VC$m3           && VC$margins[2] %in% VC$m3)           {S.h <- adiag(EQ2P, EQ4P, EQ6P)  } 
       if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% VC$m3)           {S.h <- adiag(EQ2P, EQ4P, EQ5P)  }  
       if(VC$margins[1] %in% VC$m3           && VC$margins[2] %in% c(VC$m2,VC$m2d)) {S.h <- adiag(EQ2P, EQ4P)        }     
    
    }    
    
        
        
        
}    ### TRIV    
        
        
        
        
        
        
####################################################################################################        
        
if(VC$triv == TRUE){ # here I am setting up S.h but also changing/expanding qu.mag 

S.h <- adiag(EQ1P, EQ2P, EQ3P)   

if(VC$penCor %in% c("unpen")) S.h <- adiag(S.h, matrix(0, 3, 3))        
            
            
if(VC$penCor %in% c("ridge")){

  A <- diag(c(1, 1, 1))
  if( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0) qu.mag$Ss[[1]]                   <- A
  if( VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0) qu.mag$Ss[[length(qu.mag$Ss)+1]] <- A

  S.h <- adiag(S.h, sp[length(sp)]*qu.mag$Ss[[length(qu.mag$Ss)]]) 
  
                             }
  
                                          
        
} 

####################################################################################################        

        
list(S.h = S.h, qu.mag = qu.mag)
    
   
 
         }


