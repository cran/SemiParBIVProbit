plot.SemiParBIVProbit <- function(x, eq, select, rug=TRUE, se=TRUE, se.l=1.95996, seWithMean=FALSE, n=100,
                                 xlab = NULL, ylab=NULL, zlab=NULL, xlim=NULL, ylim = NULL, main=NULL, trans = I, n2 = 40, 
                                 theta = 30, phi = 30, too.far = 0.1, ...){
  sub.edf <- function(lab, edf){ 
      pos <- regexpr(":", lab)[1]
      if(pos < 0){pos <- nchar(lab) - 1
                  lab <- paste(substr(lab, start = 1, stop = pos),",", round(edf, digits = 2), ")", sep = "")
      }
      else{lab1 <- substr(lab, start = 1, stop = pos - 2)
           lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
           lab <- paste(lab1, ",", round(edf, digits = 2), lab2,sep = "")
      }
  lab
  }

  if(x$l.sp1==0 && x$l.sp2==0) stop("The model is fully parametric; no smooth components to plot")
  
  if(eq==1) if(select>x$l.sp1) stop("No more smooth component to plot")
  if(eq==2) if(select>x$l.sp2) stop("No more smooth component to plot")
  
  if(eq==1) if(x$gam1$smooth[[select]]$dim>2) stop("No plotting for smooths of more than two variables")
  if(eq==2) if(x$gam2$smooth[[select]]$dim>2) stop("No plotting for smooths of more than two variables")

  Kk1 <- Kk2 <- 0; if(x$npRE==TRUE){ Kk1 <- x$K - 1; Kk2 <- Kk1 + 1}  

  if(eq==1) ind <- (x$gam1$smooth[[select]]$first.para+Kk1):(x$gam1$smooth[[select]]$last.para+Kk1)
  if(eq==2) ind <- (x$gam2$smooth[[select]]$first.para:x$gam2$smooth[[select]]$last.para+Kk1)+x$X1.d2+Kk2
  if(eq==1) gam <- x$gam1 else gam <- x$gam2

  est.par.c1 <- x$fit$argument
  est.par.c2 <- est.par.c1[-c(1:(x$X1.d2+Kk2))]
  
  lf  <- length(x$fit$argument)  
  Vb  <- x$Vb[1:lf,1:lf]
  d.F <- diag(x$F[1:lf,1:lf])

         edf <- round(sum(d.F[ind]),2)
         if(gam$smooth[[select]]$dim==1){raw <- gam$model[gam$smooth[[select]]$term] 
                                                 xx  <- seq(min(raw), max(raw), length=n) 
                                                      if(gam$smooth[[select]]$by!= "NA"){by <- rep(1, n)
                                                                                                 d  <- data.frame(x = xx, by = by)
                                                                                                 names(d) <- c(gam$smooth[[select]]$term,gam$smooth[[select]]$by) 
                                                      }
                                                      else{d <- data.frame(x = xx)
                                                           names(d) <- gam$smooth[[select]]$term
                                                      }
         }
         else if(gam$smooth[[select]]$dim==2){xterm <- gam$smooth[[select]]$term[1]
                                                      yterm <- gam$smooth[[select]]$term[2]      
                 						      raw <- data.frame(x = as.numeric(gam$model[xterm][[1]]),
                              			                        y = as.numeric(gam$model[yterm][[1]]) )
               						      n2 <- max(10, n2)
               						      xm <- seq(min(raw$x), max(raw$x), length = n2)
               						  	ym <- seq(min(raw$y), max(raw$y), length = n2)
              						 	xx <- rep(xm, n2)
              						      yy <- rep(ym, rep(n2, n2))
                					       	   if(too.far > 0) exclude <- exclude.too.far(xx, yy, raw$x, raw$y,dist = too.far)
                                                	   else exclude <- rep(FALSE, n2 * n2)
                                                 	   if(gam$smooth[[select]]$by != "NA"){by <- rep(1, n2^2)
                                                 	                                               d  <- data.frame(x = xx, y = yy, by = by)
                                                 	                                               names(d) <- c(xterm, yterm, gam$smooth[[select]]$by)
                                                  	   }
                                                   	   else{d <- data.frame(x = xx, y = yy)
                                                    	        names(d) <- c(xterm, yterm)
                                                         } 
              }
  X <- PredictMat(gam$smooth[[select]], d) 
  if(eq==1) f <- X%*%est.par.c1[ind] else f <- X%*%est.par.c2[ind-(x$X1.d2+Kk2)]
  

  if(se){
      if(seWithMean==TRUE && x$npRE==FALSE){
                if(eq==1) X1 <- matrix(c(gam$cmX,rep(0,length(est.par.c2))), nrow(X), ncol(Vb), byrow = TRUE) else X1 <- matrix(c(rep(0,length(gam$cmX)),gam$cmX,0), nrow(X), ncol(Vb), byrow = TRUE) 
                X1[,ind] <- X
                se.fit <- sqrt(rowSums((X1 %*% Vb) * X1))                           
                                           }
      else se.fit <- sqrt(rowSums((X %*% Vb[ind,ind]) * X))
         
  ub <- (f+se.l*se.fit)
  lb <- (f-se.l*se.fit) 
  if(eq==1 && gam$smooth[[select]]$dim==1) if(is.null(ylim)) ylim <- c(min(lb),max(ub))
  if(eq==2 && gam$smooth[[select]]$dim==1) if(is.null(ylim)) ylim <- c(min(lb),max(ub))

  } 
              


         if(gam$smooth[[select]]$dim==1){
       	   if(is.null(xlab)) x.lab <- gam$smooth[[select]]$term else x.lab <- xlab  
       	   if(is.null(ylab)) y.lab <- sub.edf(gam$smooth[[select]]$label,edf) else y.lab <- ylab  

      	 plot(xx,f,type="l",xlab=x.lab,ylim=ylim,ylab=y.lab, main=main, ...)
      	   if(se){lines( xx , ub, lty=2 )
       	          lines( xx , lb, lty=2 )
       	   }
      	   if(rug) rug(as.numeric(gam$model[gam$smooth[[select]]$term][[1]]))
         }
         else{if(gam$smooth[[select]]$dim==2){     if(is.null(zlab)) zlabel <- sub.edf(gam$smooth[[select]]$label, edf) else zlabel <- zlab
                                                      if(is.null(xlab)) xlabel <- xterm else xlabel <- xlab
                                                      if(is.null(ylab)) ylabel <- yterm else ylabel <- ylab


                pd <- list(fit = f, dim = 2, xm = xm, ym = ym, ylab = ylabel, xlab = xlabel, zlab = zlabel, raw = raw)
                if (is.null(ylim)) pd$ylim <- range(ym) else pd$ylim <- ylim
                if (is.null(xlim)) pd$xlim <- range(xm) else pd$xlim <- xlim
                         persp(pd$xm, pd$ym, matrix(trans(pd$fit), n2, n2), xlab = pd$xlab, 
                               ylab = pd$ylab, zlab = pd$zlab, 
                               ylim = pd$ylim, xlim = pd$xlim,
                               theta = theta, phi = phi, main = main, ...)
                
              }
          }
 
          

}
