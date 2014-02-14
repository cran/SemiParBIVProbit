S.mPL <- function(qu.mag,eqPL,start.v){

  	if( !(is.null(qu.mag)) ){


  		if(eqPL=="both"){ qu.mag$off  <- c(qu.mag$off, which( names(start.v) %in% c("xi1.star","xi2.star"))-1 )
                 		  qu.mag$rank <- c(qu.mag$rank,1,1) 
                                  lS <- length(qu.mag$Ss)
                                  qu.mag$Ss[[lS+1]] <- qu.mag$Ss[[lS+2]] <- matrix(1,1,1)  
                                }

                if(eqPL=="first"){ qu.mag$off  <- c(qu.mag$off,  which( names(start.v) %in% c("xi1.star"))-1  )
                                   qu.mag$rank <- c(qu.mag$rank,1) 
                                   lS <- length(qu.mag$Ss)
                                   qu.mag$Ss[[lS+1]] <- matrix(1,1,1)  
                                 }
                if(eqPL=="second"){ qu.mag$off  <- c(qu.mag$off, which( names(start.v) %in% c("xi2.star"))-1 )
                                    qu.mag$rank <- c(qu.mag$rank,1) 
                                    lS <- length(qu.mag$Ss)
                                    qu.mag$Ss[[lS+1]] <- matrix(1,1,1)  
                                  }
        }else{
  
        qu.mag$Ss <- list()

  		if(eqPL=="both"){ qu.mag$off  <- which( names(start.v) %in% c("xi1.star","xi2.star"))-1
                 		  qu.mag$rank <- c(1,1) 
                                  qu.mag$Ss[[1]] <- qu.mag$Ss[[2]] <- matrix(1,1,1)  
                                }

                if(eqPL=="first"){ qu.mag$off  <- which( names(start.v) %in% c("xi1.star"))-1
                                   qu.mag$rank <- 1 
                                   qu.mag$Ss[[1]] <- matrix(1,1,1)  
                                 }
                if(eqPL=="second"){ qu.mag$off  <- which( names(start.v) %in% c("xi2.star"))-1
                                    qu.mag$rank <- 1 
                                    qu.mag$Ss[[1]] <- matrix(1,1,1)  
                                  }

             }


qu.mag


}


