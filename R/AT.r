AT <- function(x, nm.end, E = TRUE, treat = TRUE, type = "bivariate", ind = NULL, delta = FALSE, 
   n.sim = 100, prob.lev = 0.05, hd.plot = FALSE, prob.plot = FALSE, 
   main = "Histogram and Kernel Density of Simulated Average Effects", 
   xlab = "Simulated Average Effects", ...){

etap.noi <- X.int <- X.noi <- eti1 <- eti0 <- etno <- indS <- bs <- ind.excl <- p.int1 <- p.int0 <- d.int1 <- d.int0 <- p.etn <- d.etn <- ass.p <- ass.pst <- C.11 <- C.10 <- sig2 <- peti1s <- peti0s <- sigma2.st <- sigma2s <- eti1s <- eti0s <- d0 <- d1 <- p.etns <- etnos <- etds <- ass.ps <- 1
diffEf <- fy1.y2 <- est.ATso <- y2 <- CIF <- Pr <- Effects <- C.11 <- C.10 <- NULL

m2  <- c("N","GU","rGU","LO","LN","WEI","WEI2","iG","GA")
m3  <- c("DAGUM")
mat <- c("DAGUM","GU","rGU","LO","LN","WEI","WEI2","iG","GA")

end <- 0
epsilon <- 0.0000001 # 0.9999999 sqrt(.Machine$double.eps)
max.p   <- 0.9999999
est.ATb <- NA
indD <- list()

if(x$ig[[1]]$response %in% x$ig[[2]]$pred.names ) {end <- 1; eq <- 2} 
if(x$ig[[2]]$response %in% x$ig[[1]]$pred.names ) {end <- 2; eq <- 1}    

if(x$margins[2] == "DAGUM" && eq == 2) { if( min(sqrt(x$sigma2)) <= 1) stop("Sigma parameter has value(s) smaller than 1, hence the mean is indeterminate.")}


if( !( type %in% c("naive","univariate","bivariate") ) ) stop("Error in parameter type value. It should be one of: naive, univariate or bivariate.")
if(x$margins[2] != "probit" && delta == TRUE) stop("Delta method option not available.")
if(missing(nm.end)) stop("You must provide the name of the endogenous variable.")

if(x$Model=="BSS" || x$Model=="BPO" || x$Model=="BPO0" || end==0) stop("Calculation of this average treatment effect is valid for recursive models only.")

if(type == "univariate" && x$margins[2] == "N" && eq == 2 && x$gamlssfit == FALSE) stop("You need to fit the univariate model to obtain the AT. Refit the model and set gamlssfit = TRUE.")

if(x$margins[2] %in% mat && eq == 2) stop("AT currently available for Gaussian outcome margin only.")

if(is.character(nm.end)==FALSE) stop("nm.end is not a character!")
if( !is.null(ind) && E == FALSE) stop("ind is not designed to be used when some observations are excluded from the AT's calculation.")
if( type == "naive" && E == FALSE) stop("It does not make sense to calculate the naive estimate from the treated only.")

if( !is.null(ind) ){ 

    if(is.logical(ind) == FALSE) stop("ind must be a logical variable.")
    if(length(table(ind))!=2 )   stop("ind must be a logical binary variable.")
    if( length(ind) != x$n )     stop("ind must have the same length as the number of observations used in fitting.")   

}


######################################################################


if(is.null(ind)) ind <- 1:x$n


if(E == FALSE ) {

 if(x$margins[2] != "probit") ind <- 1:x$n  
 
 if(eq==1) X.int <- as.matrix(x$X1[ind,])
 if(eq==2) X.int <- as.matrix(x$X2[ind,]) 

    if(treat == TRUE)  ind <- as.logical(X.int[, nm.end]) 
    if(treat == FALSE) ind <- as.logical(X.int[, nm.end])!=TRUE
                                              
}




######################################################################

if(type == "naive" && x$margins[2] != "probit") stop("Please fit a bivariate model with intercept and endogenous variable only and then use AT with the univariate type option.")

######################################################################
######################################################################

if(type == "naive" && x$margins[2] == "probit"){ ## binary binary case with eq = 1 or eq = 2

if(eq==2){
y1 <- x$y1[ind] 
y2 <- x$y2[ind]
}

if(eq==1){
y1 <- x$y2[ind] 
y2 <- x$y1[ind]
}

tab2 <- table(y1, y2)

pY1cT1 <- prop.table(tab2,1)[4] 
pY1cT0 <- prop.table(tab2,1)[3] 

est.AT <- (pY1cT1 - pY1cT0)

sv <- qnorm(prob.lev/2,lower.tail = FALSE) * sqrt( (pY1cT1*(1-pY1cT1))/x$n + (pY1cT0*(1-pY1cT0))/x$n )

CIs <- c(est.AT - sv, est.AT + sv)

est.ATb <- est.ATso <- NULL

}

######################################################################
######################################################################

if(type != "naive" && x$margins[2] == "probit"){ ## binary binary case with eq = 1 or eq = 2

#############

if(type == "bivariate"){
	indD[[1]] <- 1:x$X1.d2 
	indD[[2]] <- x$X1.d2+(1:x$X2.d2)
                       }

if(eq==1){ X.int <- as.matrix(x$X1[ind,])
    if(type == "bivariate") ind.int <- indD[[1]]
         }

if(eq==2){ X.int <- as.matrix(x$X2[ind,])
    if(type == "bivariate") ind.int <- indD[[2]] 
         }

if(type == "bivariate") coef.int <- as.numeric(coef(x)[ind.int])
	   
             
d0 <- d1 <- X.int
d0[,nm.end] <- 0
d1[,nm.end] <- 1


if(type == "bivariate"){
	eti1 <- d1%*%coef.int 
	eti0 <- d0%*%coef.int 
                       }

if(type == "univariate"){
	if(eq==1) ngam <- x$gam1
	if(eq==2) ngam <- x$gam2

	eti1 <- d1%*%coef(ngam) 
	eti0 <- d0%*%coef(ngam) 
                         }

#############

p.int1 <- pmax(pnorm(eti1), epsilon ) 
p.int1 <- ifelse(p.int1 > max.p,max.p,p.int1) 
p.int0 <- pmax(pnorm(eti0), epsilon ) 
p.int0 <- ifelse(p.int0 > max.p,max.p,p.int0) 

est.ATso <- (p.int1 - p.int0) 

est.AT <- mean( est.ATso, na.rm = TRUE )

#############

if(delta == TRUE){

   if(type == "univariate") var <- ngam$Vp
   if(type == "bivariate")  var <- x$Vb[ind.int,ind.int]
   
   dATT <- colMeans(  ( c(dnorm(eti1))*d1 - c(dnorm(eti0))*d0 )  )
   delta.AT <- sqrt( t(dATT)%*%var%*%dATT )
                     
                 }

if(delta == FALSE){

 if(type == "univariate") {bs <- rMVN(n.sim, mean = coef(ngam), sigma=ngam$Vp); eti1s <- d1%*%t(bs);           eti0s <- d0%*%t(bs) }
 if(type == "bivariate")  {bs <- rMVN(n.sim, mean = coef(x), sigma=x$Vb);       eti1s <- d1%*%t(bs[,ind.int]); eti0s <- d0%*%t(bs[,ind.int]) } 

 peti1s <- pmax(pnorm(eti1s), epsilon )  
 peti1s <- ifelse(peti1s > max.p,max.p,peti1s)  
 peti0s <- pmax(pnorm(eti0s), epsilon ) 
 peti0s <- ifelse(peti0s > max.p,max.p,peti0s)  
 est.ATb <- colMeans(  (peti1s - peti0s) , na.rm = TRUE    ) 

                   }


}

######################################################################
######################################################################


if(type != "naive" && x$margins[2] != "probit" && eq == 1){


 y2   <- c( min(floor(x$y2)) : max(ceiling(x$y2)) ) 
 ly2  <- length(y2)
 data <- x$dataset[ind,] 
 
 fy1.y2 <- est.ATb <- NULL
 nu <- 1
 fy1.y2S <- matrix(NA, n.sim, ly2)
 diffEfS <- matrix(NA, n.sim, ly2 - 1) 

 ind.int <- 1:x$X1.d2 
 
 if(type == "bivariate"){
              bs <- rMVN(n.sim, mean = coef(x), sigma = x$Vb)
              coefe  <- x$coef[ind.int]
              coefes <- t(bs[, ind.int])
 
 }
 
 if(type == "univariate"){ 
              bs <- rMVN(n.sim, mean = coef(x$gam1), sigma = x$gam1$Vp)
              coefe <- x$gam1$coefficient
              coefes <- t(bs)
 }             
 

###########

for(i in 1:ly2) { 

	data[, 2] <- y2[i]

	eta1 <- as.matrix( predict.gam(x$gam1, newdata = data, type = "lpmatrix") )%*%coefe  

	p1 <- pnorm(eta1) 
	fy1.y2[i] <- mean(p1)

               }

                                                  
 for(j in 1:ly2){
 
   data[, 2] <- y2[j] 
   
   etins <- as.matrix( predict.gam(x$gam1, newdata = data, type = "lpmatrix") )%*%coefes
                  
   p1s <- pnorm(etins)    
   fy1.y2S[, j] <- apply(p1s, MARGIN = 2, FUN = mean, na.rm=TRUE)
                    
                 }

###########

diffEf <- diff(fy1.y2)
est.AT <- mean(diffEf)

diffEfS <- t(apply(fy1.y2S, MARGIN = 1, FUN = diff)) 
est.ATb <- apply(diffEfS, MARGIN = 1, FUN = mean) 

CIF    <- t(apply(fy1.y2S, MARGIN=2, FUN=quantile, probs=c(prob.lev/2,1-prob.lev/2), na.rm=TRUE ))
CIDiff <- t(apply(diffEfS, MARGIN=2, FUN=quantile, probs=c(prob.lev/2,1-prob.lev/2), na.rm=TRUE ))

Pr <- data.frame(Pr = fy1.y2, CIF)
names(Pr)[2:3] <- dimnames(CIF)[[2]]
dimnames(Pr)[[1]] <- y2

Effects <- data.frame(Effects = diffEf, CIDiff)  
names(Effects)[2:3] <- dimnames(CIDiff)[[2]]
dimnames(Effects)[[1]] <- y2[2:ly2]


if(prob.plot == TRUE){

plot(y2, fy1.y2, ylab = "Pr(Response=1|Treatment)", xlab = "Treatment", pch = 16, ylim = c(min(CIF[,1]),max(CIF[,2])))
lines(y2, fy1.y2, type = "l")
for (i in 1:ly2) lines( y = c(CIF[i,1], CIF[i,2]), x = c(y2[i],y2[i]) )

if(hd.plot == TRUE) dev.new()

                     }


}






if(type != "naive" && x$margins[2] == "N" && eq == 2){

 if(type == "univariate") {bs <- rMVN(n.sim, mean = x$gamlss$fit$argument, sigma=x$gamlss$magpp$Vb)
                           est.AT  <- est.ATso <- x$gamlss$fit$argument[nm.end] 
                           est.ATb <- bs[, which(names(x$gamlss$fit$argument)==nm.end) ]
                           } 
                           
 if(type == "bivariate")  {bs <- rMVN(n.sim, mean = coef(x), sigma=x$Vb)
                           est.AT  <- est.ATso <- x$coefficients[nm.end]
                           est.ATb <- bs[, nm.end]        
                           }

}
  


######################################################################
######################################################################


  if(hd.plot == TRUE && type != "naive"){
  
  if(x$margins[2] %in% c(m2, m3) && eq == 2) mult <- 1 else mult <- 100
  
  hist(est.ATb*mult, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(est.ATb*mult)$y,hist(est.ATb*mult, plot = FALSE)$density)), ...)
  lines(density(est.ATb*mult))

                                         }



if( !(type == "naive" && x$margins[2] == "probit") ){

if(type != "naive" && x$margins[2] == "probit" && delta == TRUE){ 

                   esd.AT <- delta.AT*qnorm(prob.lev/2,lower.tail = FALSE) 
                   CIs <- c(est.AT - esd.AT, est.AT + esd.AT)
                   
               }else CIs <- as.numeric(quantile(est.ATb, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))


}







res <- c(CIs[1], est.AT, CIs[2])

out <- list(res=res, prob.lev=prob.lev, sim.AT=est.ATb, AT.so = est.ATso, mar2=x$margins[2], type = type, 
            Effects = Effects, Pr = Pr, treat = y2, eq = eq)
 							 
rm(etap.noi, X.int, X.noi, eti1, eti0, etno, indS, bs, ind.excl, p.int1, p.int0, d.int1, d.int0,
   p.etn, d.etn, ass.p, ass.pst, C.11, C.10, sig2, peti1s, peti0s, sigma2.st, sigma2s, eti1s, eti0s, d0, d1,
   p.etns, etnos, etds, ass.ps)    
 
class(out) <- "AT"

out





}





















