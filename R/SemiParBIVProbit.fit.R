SemiParBIVProbit.fit <- function(func.opt, start.v, rinit, rmax, BivD, nC, nu, PL, eqPL, H.n, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, control.sp, gamma, X1.d2, X2.d2, pPen1=NULL, pPen2=NULL, sp=NULL, qu.mag=NULL, gp1, gp2, fp, aut.sp, l.sp1, l.sp2, pr.tolsp, weights, iterlimsp, 
fterm, mterm, iterlim){ 

spE <- sp

if( PL!="P" && (l.sp1!=0 || l.sp2!=0) ){
    lspe <- length(sp)
    if(eqPL=="both") exclu <- c(lspe-1,lspe) else exclu <- c(lspe)  
    spE <- sp[-exclu]
    qu.mag$exclu <- exclu
}

  fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, BivD=BivD, nC=nC, nu=nu, sp.xi1=sp["xi1"], sp.xi2=sp["xi2"], PL=PL, eqPL=eqPL, H.n=H.n, 
                y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                X1=X1, X2=X2,  
                X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=spE, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, blather=TRUE, weights=weights, 
                fterm=fterm, mterm=mterm, iterlim = iterlim)  

  iter.if <- fit$iterations  

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- NULL

       
    if((aut.sp==TRUE && fp==FALSE && (l.sp1!=0 || l.sp2!=0)) || PL!="P"){

       stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0  
       
       
         if(PL == "P")          myf <- function(x) list( rbind( c(x[1],x[2]), 
                                                                c(x[2],x[3])  ) )
       
         if(PL != "P"){
       
              if(eqPL!="both")  myf <- function(x) list( rbind( c(x[1],x[2],x[4]),
                                                                c(x[2],x[3],x[5]),
                                                                c(x[4],x[5],x[6])  ) )
                    
              if(eqPL=="both")  myf <- function(x) list( rbind( c(x[1],x[2],x[4],x[5]),
                                                                c(x[2],x[3],x[6],x[7]),
                                                                c(x[4],x[6],x[8],x[9]),
                                                                c(x[5],x[7],x[9],x[10])  ) )
                      }
       


	  while( stoprule.SP > pr.tolsp ){ 

             coefo <- fit$argument
             spEo <- spo <- sp 
             if(PL!="P" && (l.sp1!=0 || l.sp2!=0) ) spEo <- spo[-exclu]   
  
		 wor.c <- try(working.comp(fit,X1,X2,X1.d2,X2.d2,myf)) 
                 if(class(wor.c)=="try-error") break
             
                	bs.mgfit <- try(magic(y=wor.c$rW.Z,X=wor.c$rW.X,sp=sp,S = qu.mag$Ss,
                        	           off=qu.mag$off,rank=qu.mag$rank,
                                	   gcv=FALSE,gamma=gamma,control=control.sp))
                		if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                	spE <- sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1; names(sp) <- names(spo) 
                        if(PL!="P" && (l.sp1!=0 || l.sp2!=0) ) spE <- sp[-exclu]
                           
             o.ests <- c(fit$argument)

             fit <- try(trust(func.opt, fit$argument, rinit=rinit, rmax=rmax, BivD=BivD, nC=nC, nu=nu, sp.xi1=sp["xi1"], sp.xi2=sp["xi2"], PL=PL, eqPL=eqPL, H.n=H.n, 
                              y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,   
                             X1=X1, X2=X2,  
                              X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=spE, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, weights=weights, blather=TRUE, 
                              iterlim=1e+4, fterm=fterm, mterm=mterm),silent=TRUE)

              			if(class(fit)=="try-error"){ 
              	 			fit  <- trust(func.opt, coefo, rinit=rinit, rmax=rmax, BivD=BivD, nu=nu, sp.xi1=spo["xi1"], sp.xi2=spo["xi2"], nC=nC, PL=PL, eqPL=eqPL, H.n=H.n,  
                                                      y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                                                      X1=X1, X2=X2,  
                		              	      X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=spEo, qu.mag=qu.mag, gp1=gp1, gp2=gp2, 
                                                      fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, blather=TRUE, weights=weights,
                        		              iterlim=1e+4, fterm=fterm, mterm=mterm)
               		        conv.sp <- FALSE; break
                	                                    } 
             iter.inner <- iter.inner + fit$iterations   	                                    
              

             n.ests <- c(fit$argument)

             if(iter.sp>iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- max(abs(o.ests-n.ests))
                  
           
          }
      
    }


                  list(fit=fit, iter.if=iter.if, conv.sp=conv.sp, iter.sp=iter.sp, 
                       iter.inner=iter.inner, bs.mgfit=bs.mgfit, wor.c=wor.c, sp=sp)


}




