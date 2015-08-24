S.m <- function(gam1, gam2, gam3, gam4, gam5, gam6, l.sp1, l.sp2, l.sp3, l.sp4, l.sp5, l.sp6, eq1 = "yes", margin2){

  m2   <- c("N","GU","rGU","LO","LN","WEI","WEI2","iG","GA")
  m3   <- c("DAGUM")

  Ss <- list()
  off <- rank <- 0 
  i1 <- i2 <- i3 <- i4 <- i5 <- i6 <- 1
  
  l.gam1 <- length(coef(gam1)) 
  
  if(eq1 == "no"){  l.gam1 <- 0 
  
    if(margin2 %in% m2)  l.sp1 <- l.sp4 <- l.sp5 <- l.sp6 <- 0 
    if(margin2 %in% m3)  l.sp1 <- l.sp5 <- l.sp6 <- 0 
 
  } 
  
	for( j in 1:(l.sp1 + l.sp2 + l.sp3 + l.sp4 + l.sp5 + l.sp6 ) ){
	
		if(j <= l.sp1){                          
		             Ss[[j]] <- gam1$smooth[[i1]]$S[[1]]
            		     rank[j] <- gam1$smooth[[i1]]$rank
                              off[j] <- gam1$smooth[[i1]]$first.para 
                                  i1 <- i1 + 1                            
                              }   
                              
                if( (j >  l.sp1 &&  j <= (l.sp1 + l.sp2) ) ){        
                     	     Ss[[j]] <- gam2$smooth[[i2]]$S[[1]]  
                             off[j]  <- l.gam1 + gam2$smooth[[i2]]$first.para 
                             rank[j] <- gam2$smooth[[i2]]$rank
                                  i2 <- i2 + 1                              
                              }     
                              
               if(j >  (l.sp1 + l.sp2) && j <= (l.sp1 + l.sp2 + l.sp3)  ){        
                     	     Ss[[j]] <- gam3$smooth[[i3]]$S[[1]]  
                             off[j]  <- l.gam1 + length(coef(gam2)) + gam3$smooth[[i3]]$first.para 
                             rank[j] <- gam3$smooth[[i3]]$rank
                                  i3 <- i3 + 1                              
                              } 
                 
                if(j >  (l.sp1 + l.sp2 + l.sp3) && j <= (l.sp1 + l.sp2 + l.sp3 + l.sp4)  ){        
                     	     Ss[[j]] <- gam4$smooth[[i4]]$S[[1]]  
                             off[j]  <- l.gam1 + length(coef(gam2)) + length(coef(gam3)) + gam4$smooth[[i4]]$first.para 
                             rank[j] <- gam4$smooth[[i4]]$rank
                                  i4 <- i4 + 1                              
                              }                  
                 
      
                if(j >  (l.sp1 + l.sp2 + l.sp3 + l.sp4) && j <= (l.sp1 + l.sp2 + l.sp3 + l.sp4 + l.sp5)  ){        
                     	     Ss[[j]] <- gam5$smooth[[i5]]$S[[1]]  
                             off[j]  <- l.gam1 + length(coef(gam2)) + length(coef(gam3)) + length(coef(gam4)) + gam5$smooth[[i5]]$first.para 
                             rank[j] <- gam5$smooth[[i5]]$rank
                                  i5 <- i5 + 1                              
                              }                            
                              
                if(j >  (l.sp1 + l.sp2 + l.sp3 + l.sp4 + l.sp5) && j <= (l.sp1 + l.sp2 + l.sp3 + l.sp4 + l.sp5 + l.sp6)  ){        
                     	     Ss[[j]] <- gam6$smooth[[i6]]$S[[1]]  
                             off[j]  <- l.gam1 + length(coef(gam2)) + length(coef(gam3)) + length(coef(gam4)) + length(coef(gam5)) + gam6$smooth[[i6]]$first.para 
                             rank[j] <- gam6$smooth[[i6]]$rank
                                  i6 <- i6 + 1                              
                              }                                
                                                           
}


       
list(rank=rank,off=off,Ss=Ss)

}


