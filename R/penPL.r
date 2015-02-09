penPL <- function(qu.mag, sp, VC, add.z, fitPL){
              
        dimP1 <- dimP2 <- 0     
        S1 <- S2 <- matrix(0,1,1)  
    
        S <- mapply("*", qu.mag$Ss[-qu.mag$exclu], sp, SIMPLIFY=FALSE)
        S <- do.call(adiag, lapply(S, unlist))
    
        ma1 <- matrix(0,VC$gp1,VC$gp1) 
        ma2 <- matrix(0,VC$gp2,VC$gp2)
    
        if(0!=0){ indP1 <- qu.mag$off[1]:(qu.mag$off[1]+qu.mag$rank[1]-1)
                              dimP1 <- length(indP1)
                              ma1[indP1,indP1] <- S[1:dimP1,1:dimP1]
                                    } 
    
        if(0!=0){ 
                              indP2 <- (qu.mag$off[VC$l.sp1+1]-VC$X1.d2):(-VC$X1.d2+qu.mag$off[VC$l.sp1+1]+qu.mag$rank[VC$l.sp1+1]-1)
                              dimP2 <- length(indP2)
                              ma2[indP2,indP2] <- S[(dimP1+1):(length(indP2)+dimP1),(dimP1+1):(length(indP2)+dimP1)]
                                    }                                 
        
        lP1 <- 0; lP2 <- 0
        
        if((lP1!=0 && VC$l.sp1>1) || (lP1==0 && VC$l.sp1>0)) S1 <- S[(dimP1+1):(dimP1+VC$X1.d2-VC$gp1),(dimP1+1):(dimP1+VC$X1.d2-VC$gp1)]
        if((lP2!=0 && VC$l.sp2>1) || (lP2==0 && VC$l.sp2>0)){dS1 <- dim(S1)[2]; if(dS1==1) dS1 <- 0; 
                                                       S2 <- S[(dimP1+dimP2+dS1+1):dim(S)[2],(dimP1+dimP2+dS1+1):dim(S)[2]]}
        
        lS1 <- length(S1); lS2 <- length(S2) 
    
    if(fitPL!="fixed"){
    
        if(lS1==1 && lS2==1) S.h <- adiag(ma1, ma2, 0, add.z)
        if(lS1 >1 && lS2==1) S.h <- adiag(ma1, S1, ma2, 0, add.z)
        if(lS1==1 && lS2 >1) S.h <- adiag(ma1, ma2, S2, 0, add.z)
        if(lS1 >1 && lS2 >1) S.h <- adiag(ma1, S1, ma2, S2, 0, add.z)
            
      }else{
      
        if(lS1==1 && lS2==1) S.h <- adiag(ma1, ma2, 0)
        if(lS1 >1 && lS2==1) S.h <- adiag(ma1, S1, ma2, 0)
        if(lS1==1 && lS2 >1) S.h <- adiag(ma1, ma2, S2, 0)
        if(lS1 >1 && lS2 >1) S.h <- adiag(ma1, S1, ma2, S2, 0)
      
      }
            
            
             
             
list(S.h = S.h)
   
         }





















