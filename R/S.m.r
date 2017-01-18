S.m <- function(GAM, L.SP, L.GAM){

  Ss <- list()
  off <- rank <- 0 
  i1 <- i2 <- i3 <- i4 <- i5 <- i6 <- i7 <- i8 <- 1
  
	for( j in 1:(L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6 + L.SP$l.sp7 + L.SP$l.sp8) ){
	
		if(j <= L.SP$l.sp1){                          
		             Ss[[j]] <- GAM$gam1$smooth[[i1]]$S[[1]]
            		     rank[j] <- GAM$gam1$smooth[[i1]]$rank
                              off[j] <- GAM$gam1$smooth[[i1]]$first.para 
                                  i1 <- i1 + 1                            
                              }   
                              
                if( (j >  L.SP$l.sp1 &&  j <= (L.SP$l.sp1 + L.SP$l.sp2) ) ){        
                     	     Ss[[j]] <- GAM$gam2$smooth[[i2]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + GAM$gam2$smooth[[i2]]$first.para 
                             rank[j] <- GAM$gam2$smooth[[i2]]$rank
                                  i2 <- i2 + 1                              
                              }     
                              
               if(j >  (L.SP$l.sp1 + L.SP$l.sp2) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3)  ){        
                     	     Ss[[j]] <- GAM$gam3$smooth[[i3]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + GAM$gam3$smooth[[i3]]$first.para 
                             rank[j] <- GAM$gam3$smooth[[i3]]$rank
                                  i3 <- i3 + 1                              
                              } 
                 
                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4)  ){        
                     	     Ss[[j]] <- GAM$gam4$smooth[[i4]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + GAM$gam4$smooth[[i4]]$first.para 
                             rank[j] <- GAM$gam4$smooth[[i4]]$rank
                                  i4 <- i4 + 1                              
                              }                  
                 
      
                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5)  ){        
                     	     Ss[[j]] <- GAM$gam5$smooth[[i5]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + GAM$gam5$smooth[[i5]]$first.para 
                             rank[j] <- GAM$gam5$smooth[[i5]]$rank
                                  i5 <- i5 + 1                              
                              }                            
                              
                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6)  ){        
                     	     Ss[[j]] <- GAM$gam6$smooth[[i6]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + GAM$gam6$smooth[[i6]]$first.para 
                             rank[j] <- GAM$gam6$smooth[[i6]]$rank
                                  i6 <- i6 + 1                              
                              }  

                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6 + L.SP$l.sp7)  ){        
                     	     Ss[[j]] <- GAM$gam7$smooth[[i7]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + L.GAM$l.gam6 + GAM$gam7$smooth[[i7]]$first.para 
                             rank[j] <- GAM$gam7$smooth[[i7]]$rank
                                  i7 <- i7 + 1                              
                              }  
                              
                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6 + L.SP$l.sp7) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6 + L.SP$l.sp7 + L.SP$l.sp8)  ){        
                     	     Ss[[j]] <- GAM$gam8$smooth[[i8]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + L.GAM$l.gam6 + L.GAM$l.gam7 + GAM$gam8$smooth[[i8]]$first.para 
                             rank[j] <- GAM$gam8$smooth[[i8]]$rank
                                  i8 <- i8 + 1                              
                              }                               
                                                           
}


       
list(rank=rank, off=off, Ss=Ss)

}


