extraiterNP <- function(paramNP, y1, y2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, X1.d2, X2.d2, pPen1, pPen2, sp, qu.mag=NULL, gp1, gp2, fp, l.sp1, l.sp2, weights, 
                         K, n, N, cuid, uidf, masses, pr.tolsp, NGQ=NULL, dat1all=NULL, dat2all=NULL, W=NULL){
          
        nsv <- names(paramNP); count.npRE <- 0   
                   
        T.sv <- bprobNP(paramNP, y1=y1, y2=y2, 
                         y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                         X1=X1, X2=X2,  
                         X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, weights=weights, 
                         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses, NGQ=NGQ, dat1all=dat1all, dat2all=dat2all, W=W)
                       
      	  while( max(abs(T.sv$gradient)) > pr.tolsp*10000 && sum(as.numeric(T.sv$gradient=="NaN"))==0){
            
            	starm.v <- paramNP 
	    	paramNP <- paramNP - ginv(T.sv$hessian)%*%T.sv$gradient 
	    	masses <- T.sv$masses; count.npRE <- count.npRE + 1
	    	T.sv <- bprobNP(paramNP, y1=y1, y2=y2, 
                                y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                                     X1=X1, X2=X2,  
                	             X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, weights=weights, 
                        	     K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses, NGQ=NGQ, dat1all=dat1all, dat2all=dat2all, W=W)
	    	
	    	if(count.npRE > 5000) break }

            if(sum(as.numeric(T.sv$gradient=="NaN"))==0) paramNP <- as.vector(paramNP) else paramNP <- as.vector(starm.v) 

            names(paramNP) <- nsv  

            list(paramNP=paramNP,count.npRE=count.npRE,masses=masses) 
            

}                    
                                  
                   