S.m <- function(gam1,gam2,K,npRE){

  if(npRE==TRUE){k1 <- K-1; k2 <- 2*K-2} else {k1 <- k2 <- 0}
  Ss <- list()
  off <- rank <- 0 
  jj <- 1
	for(j in 1:(length(gam1$sp)+length(gam2$sp))){
		if(j<=length(gam1$sp)){Ss[[j]] <- gam1$smooth[[j]]$S[[1]]
            		           rank[j] <- gam1$smooth[[j]]$rank
                        		off[j] <- gam1$smooth[[j]]$first.para + k1
            }
            else{Ss[[j]] <- gam2$smooth[[jj]]$S[[1]]  
                 off[j]  <- gam1$smooth[[length(gam1$sp)]]$last.para + gam2$smooth[[jj]]$first.para + k2
                 rank[j] <- gam2$smooth[[jj]]$rank  
                 jj <- jj + 1 
            }                                            
       }
list(rank=rank,off=off,Ss=Ss)
}


