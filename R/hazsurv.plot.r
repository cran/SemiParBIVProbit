hazsurv.plot <- function(x, eq, newdata, type = "surv", intervals = TRUE, n.sim = 100, prob.lev = 0.05, 
                         shade = FALSE, ylim, ylab, xlab, ...){

pr <- h <- hs <- prs <- CIpr <- CIh <- NULL

if(missing(eq) && x$univar.gamlss == FALSE) stop("You must provide the equation number (either 1 or 2).")
if(x$univar.gamlss == TRUE) eq <- 1
if(missing(newdata))        stop("You have to provide a new data frame.")   
if(!is.data.frame(newdata)) stop("You have to provide a new data frame.")   
if( dim(newdata)[1] != 1 )  stop("Your data frame has to contain one row only.")
if(!(type %in% c("surv", "hazard"))) stop("The type argument can either be surv or hazard")
if(x$surv.flex == FALSE) stop("This function is only suitable for flexible survival models.")
if(missing(ylim)) ylim <- NULL   

#################################################################################################


if(eq == 1){
  ntv  <- as.character(x$formula[[1]][2])
  tv   <- seq(range(x$y1)[1], range(x$y1)[2], length.out = 100)
  indp <- 1:x$VC$X1.d2
  gob  <- x$gam1
}

if(eq == 2){
  ntv  <- as.character(x$formula[[2]][2])
  tv   <- seq(range(x$y2)[1], range(x$y2)[2], length.out = 100)
  indp <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)
  gob  <- x$gam2
}

ti <- data.frame(tv)
names(ti) <- ntv

newdata <- data.frame(ti, newdata)
Xpred <- predict(x, newdata, eq = eq, type = "lpmatrix")

params1 <- x$coef.t[indp]
eta1 <- Xpred%*%params1
pd <- probmS(eta1, x$VC$margins[eq])    
pr <- pd$pr 
  
bs <- rMVN(n.sim, mean = x$coef.t, sigma = x$Vb.t)
#bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  

if(!is.null(x$VC$mono.sm.pos)) mono.sm.pos <- x$VC$mono.sm.pos else mono.sm.pos <- c(x$VC$mono.sm.pos1, x$VC$mono.sm.pos2 + x$VC$X1.d2)  

bs[, mono.sm.pos] <- ifelse(bs[, mono.sm.pos] < 0, 0, bs[, mono.sm.pos]) 


#bs[, mono.sm.pos] <- exp(bs[, mono.sm.pos]) 

eta1s <- Xpred%*%t(bs[,indp])

pds <- probmS(eta1s, x$VC$margins[eq]) 
prs <- pds$pr 

#################################################################################################

if(type == "surv"){

if(intervals == TRUE) CIpr <- rowQuantiles(prs, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

if(is.null(ylim) && intervals == TRUE) ylim <- c(min(CIpr[,1]),max(CIpr[,2]) )  
if(missing(ylab))                      ylab <- "Survival function"  
if(missing(xlab))                      xlab <- "Time"  

plot(tv, pr, type = "l", ylab = ylab, xlab = xlab, ylim = ylim, ...)

if(intervals == TRUE){
  if(shade == FALSE){
    lines(tv, CIpr[,1], lty = 2)
    lines(tv, CIpr[,2], lty = 2)
                    }else{
    polygon(c(tv,rev(tv)), c(CIpr[,1],rev(CIpr[,2])), col= "gray80", border = NA)
    lines(tv, pr, type = "l")
                          }
                      }
}


#################################################################################################

if(type == "hazard"){

h <- .Machine$double.eps^(1/3) * ifelse(abs(newdata[, ntv ]) > 1, abs(newdata[, ntv ]), 1)
temp <- newdata[, ntv ] + h
h.hi <- temp - newdata[, ntv ]
temp <- newdata[, ntv ] - h
h.lo <- newdata[, ntv ] - temp
twoeps <- h.hi + h.lo
newdata1 <- newdata2 <- newdata
newdata1[, ntv ] <- newdata1[, ntv ] - h
newdata2[, ntv ] <- newdata2[, ntv ] + h
attr(newdata1, "terms") <- attr(newdata2, "terms") <- NULL 
Xd   <- (predict(gob, newdata2, type = "lpmatrix") - predict(gob, newdata1, type = "lpmatrix"))/twoeps 
rm(newdata1, newdata2)

Xthe <- Xd%*%params1   
  
Gp <- pd$dS
h  <- -Gp/pr*Xthe

Gps <- pds$dS
Xthes <- Xd%*%t(bs[,indp]) 

hs <- -Gps/prs*Xthes

if(intervals == TRUE) CIh <- rowQuantiles(hs, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

if(is.null(ylim) && intervals == TRUE) ylim <- c(min(CIh[,1]),max(CIh[,2]) )  
if(missing(ylab))                      ylab <- "Hazard"  
if(missing(xlab))                      xlab <- "Time"  

plot(tv, h, type = "l", ylab = ylab, xlab = xlab, ylim = ylim, ...)

if(intervals == TRUE){
  if(shade == FALSE){
    lines(tv, CIh[,1], lty = 2)
    lines(tv, CIh[,2], lty = 2)
                    }else{
    polygon(c(tv,rev(tv)), c(CIh[,1],rev(CIh[,2])), col= "gray80", border = NA)
    lines(tv, h, type = "l")
                         }
                     }
}

#################################################################################################

out.r <- list(s = pr, h = h, h.sim = hs, s.sim = prs)
invisible(out.r)

}

