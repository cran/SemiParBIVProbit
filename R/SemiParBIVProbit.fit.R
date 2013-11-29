SemiParBIVProbit.fit <- function(func.opt, start.v, rinit, rmax, BivD, nC, nu, xi1, xi2, PL, eqPL, H.n, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, RE, RE.type, e.npRE, control.sp, gamma,
                X1.d2, X2.d2, pPen1=NULL, pPen2=NULL, sp=NULL, qu.mag=NULL, gp1, gp2, fp, aut.sp, l.sp1, l.sp2, pr.tolsp, weights, iterlimsp, 
                fterm, mterm, iterlim, K=NULL, n=NULL, N=NULL, cuid=NULL, uidf=NULL, masses=NULL, NGQ=NULL, dat1all=NULL, dat2all=NULL, W=NULL){ 


  fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, BivD=BivD, nC=nC, nu=nu, xi1=xi1, xi2=xi2, PL=PL, eqPL=eqPL, H.n=H.n, 
                y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                X1=X1, X2=X2,  
                X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, blather=TRUE, weights=weights, 
                fterm=fterm, mterm=mterm, iterlim = iterlim, K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses, NGQ=NGQ, dat1all=dat1all, dat2all=dat2all, W=W)  

  iter.if <- fit$iterations  

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- NULL
  if(RE==FALSE || (RE==TRUE && RE.type=="N")) masses <- NULL
       
       
    if(aut.sp==TRUE && fp==FALSE && (l.sp1!=0 || l.sp2!=0) ){

       stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0  

	  while( stoprule.SP > pr.tolsp ){ 

             coefo <- fit$argument
             spo <- sp    
  
		 if(RE!=TRUE) wor.c <- try(working.comp(fit,X1,X2,X1.d2,X2.d2)) else wor.c <- try(working.compNP(fit,X1,X2,X1.d2,X2.d2,K)) 
                 if(class(wor.c)=="try-error") break
             
                	bs.mgfit <- try(magic(y=wor.c$rW.Z,X=wor.c$rW.X,sp=sp,S = qu.mag$Ss,
                        	           off=qu.mag$off,rank=qu.mag$rank,
                                	   gcv=FALSE,gamma=gamma,control=control.sp))
                		if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                	sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1 
                               
             o.ests <- c(fit$argument)
            
               if(RE==TRUE && RE.type=="NP" && e.npRE==TRUE){
		     up.ve <- extraiterNP(paramNP=o.ests, BivD, nC, nu, xi1, xi2, PL, eqPL, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, X1.d2, X2.d2, pPen1, pPen2, sp, qu.mag, gp1, gp2, fp, l.sp1, l.sp2, weights, 
                	                     K, n, N, cuid, uidf, masses, pr.tolsp, NGQ, dat1all, dat2all, W) 
                     fit$argument <- up.ve$paramNP; masses <- up.ve$masses
                             		     }

             fit <- try(trust(func.opt, fit$argument, rinit=rinit, rmax=rmax, BivD=BivD, nC=nC, nu=nu, xi1=xi1, xi2=xi2, PL=PL, eqPL=eqPL, H.n=H.n, 
                              y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,   
                             X1=X1, X2=X2,  
                              X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, weights=weights, blather=TRUE, 
                              iterlim=1e+4, fterm=fterm, mterm=mterm, K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses,NGQ=NGQ, dat1all=dat1all, dat2all=dat2all, W=W),silent=TRUE)

              			if(class(fit)=="try-error"){ 
              	 			fit  <- trust(func.opt, coefo, rinit=rinit, rmax=rmax, BivD=BivD, nu=nu, xi1=xi1, xi2=xi2, nC=nC, PL=PL, eqPL=eqPL, H.n=H.n,  
                                                      y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                                                      X1=X1, X2=X2,  
                		              	      X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=spo, qu.mag=qu.mag, gp1=gp1, gp2=gp2, 
                                                      fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, blather=TRUE, weights=weights,
                        		              iterlim=1e+4, fterm=fterm, mterm=mterm, K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses,NGQ=NGQ, dat1all=dat1all, dat2all=dat2all, W=W)
               		        conv.sp <- FALSE; break
                	                                    } 
             iter.inner <- iter.inner + fit$iterations   	                                    
              

             n.ests <- c(fit$argument)

             if(iter.sp>iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- max(abs(o.ests-n.ests))
                  
           
          }
      
    }


                  list(fit=fit, iter.if=iter.if, conv.sp=conv.sp, iter.sp=iter.sp, 
                       iter.inner=iter.inner, bs.mgfit=bs.mgfit, wor.c=wor.c, sp=sp, 
                       masses=masses)


}




