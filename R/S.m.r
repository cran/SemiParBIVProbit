S.m <- function(gam1, gam2, gam3, l.sp1, l.sp2, l.sp3){

  Ss <- list()
  off <- rank <- 0 
  i1 <- i2 <- i3 <- 1
  
	for( j in 1:(l.sp1 + l.sp2 + l.sp3) ){
	
		if(j <= l.sp1){                          
		             Ss[[j]] <- gam1$smooth[[i1]]$S[[1]]
            		     rank[j] <- gam1$smooth[[i1]]$rank
                              off[j] <- gam1$smooth[[i1]]$first.para 
                                  i1 <- i1 + 1                            
                              }   
                              
                if( (j >  l.sp1 &&  j <= (l.sp1 + l.sp2) ) ){        
                     	     Ss[[j]] <- gam2$smooth[[i2]]$S[[1]]  
                             off[j]  <- length(coef(gam1)) + gam2$smooth[[i2]]$first.para 
                             rank[j] <- gam2$smooth[[i2]]$rank
                                  i2 <- i2 + 1                              
                              }     
                              
                if(j >  (l.sp1 + l.sp2) ){        
                     	     Ss[[j]] <- gam3$smooth[[i3]]$S[[1]]  
                             off[j]  <- length(coef(gam1)) + length(coef(gam2)) + gam3$smooth[[i3]]$first.para 
                             rank[j] <- gam3$smooth[[i3]]$rank
                                  i3 <- i3 + 1                              
                              } 
                              
}
       
list(rank=rank,off=off,Ss=Ss)

}


