S.m <- function(gam1,gam2,l.sp1,l.sp2,K,RE,RE.type){

  if(RE==FALSE) k1 <- k2 <- 0
  if(RE==TRUE && RE.type=="NP") {k1 <- K-1; k2 <- 2*K-2}
  if(RE==TRUE && RE.type=="N") {k1 <- 0; k2 <- 1}

  Ss <- list()
  
  off <- rank <- 0 
  
  jj <- 1
  
	for( j in 1:(l.sp1 + l.sp2) ){
	
		if(j<=l.sp1){Ss[[j]] <- gam1$smooth[[j]]$S[[1]]
            		     rank[j] <- gam1$smooth[[j]]$rank
                              off[j] <- gam1$smooth[[j]]$first.para + k1
                             }
                else{
                     Ss[[j]] <- gam2$smooth[[jj]]$S[[1]]  
                     off[j]  <- length(coef(gam1)) + gam2$smooth[[jj]]$first.para + k2
                     rank[j] <- gam2$smooth[[jj]]$rank
                 
                jj <- jj + 1 
            }                                            
       }
       
       
       
list(rank=rank,off=off,Ss=Ss)



}


