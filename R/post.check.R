post.check <- function(x, main = "Histogram and Density Estimate of Residuals", 
                          main2 = "Histogram and Density Estimate of Residuals",
                           xlab = "Quantile Residuals", xlab2 = "Quantile Residuals", test = FALSE, ...){

if(x$triv == TRUE ) stop("This function is not suitable for trivariate probit models.")
y1m <- y2m <- NA

if(x$Cont == "NO"){


y2     <- x$y2
eta2   <- x$eta2
sigma2 <- x$sigma2
nu     <- x$nu



if(x$VC$margins[2] %in% c("ZTP")){
     
    ly2 <- length(y2)
    y2m <- list()
    my2 <- max(y2)
    for(i in 1:ly2){ y2m[[i]] <- seq(0, y2[i]); length(y2m[[i]]) <- my2+1} 
    y2m <- do.call(rbind, y2m)     
     
}





if(x$VC$ccss == "yes"){

eta2 <- eta2[x$inde]

if(!is.null(x$X3) && !(x$VC$margins[2] %in% x$VC$m1d) ){
sigma2 <- sigma2[x$inde]
nu     <- nu[x$inde]
                  }
                      }
                      
                      


if(x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3))    p  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2])$p2 

if(x$VC$margins[2] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){ 
                                               pd <- distrHsATDiscr(y2, eta2, sigma2, nu, x$margins[2], y2m = y2m) 
                                               p <- pd$p2
                                               d <- pd$pdf2  } 

if(x$VC$margins[2] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)) qr <- qnorm( runif(y2, p - d, p) )
if(x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3))            qr <- qnorm(p)

if(x$VC$ccss == "yes") qr <- qr - mean(qr) 


par(mfrow = c(1, 2))

hist(qr, freq = FALSE, #ylim=c(0, max(qr, hist(qr, plot = FALSE)$density) ),
     main=main,
     xlab=xlab, ylab = "Density", ...)
lines(density(qr, adjust = 2),lwd=2)

qqnorm(qr)
abline(0, 1, col = "red")





if(test == TRUE){

#print(ks.test(qr, "pnorm"))
print(shapiro.test(qr))

}



}




if(x$Cont == "YES"){

y1 <- x$y1
y2 <- x$y2


if(x$VC$margins[1] %in% c("ZTP")){
     
    ly1 <- length(y1)
    y1m <- list()
    my1 <- max(y1)
    for(i in 1:ly1){ y1m[[i]] <- seq(0, y1[i]); length(y1m[[i]]) <- my1+1} 
    y1m <- do.call(rbind, y1m)     
     
}


if(x$VC$margins[2] %in% c("ZTP")){
     
    ly2 <- length(y2)
    y2m <- list()
    my2 <- max(y2)
    for(i in 1:ly2){ y2m[[i]] <- seq(0, y2[i]); length(y2m[[i]]) <- my2+1} 
    y2m <- do.call(rbind, y2m)     
     
}









if(x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3))   p1 <- distrHsAT(x$y1, x$eta1, x$sigma21, x$nu1, x$margins[1])$p2 
if(x$VC$margins[1] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){

pd <- distrHsATDiscr(x$y1, x$eta1, x$sigma21, x$nu1, x$margins[1], y2m = y1m) 
p <- pd$p2
d <- pd$pdf2   

p1 <- runif(y1, p - d, p) 

}



if(x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3))   p2 <-      distrHsAT(x$y2, x$eta2, x$sigma22, x$nu2, x$margins[2])$p2 
if(x$VC$margins[2] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){

pd <- distrHsATDiscr(x$y2, x$eta2, x$sigma22, x$nu2, x$margins[2], y2m = y2m)
p <- pd$p2
d <- pd$pdf2   

p2 <- runif(y2, p - d, p) 

}


par(mfrow = c(2, 2))

qr1 <- qnorm(p1)
hist(qr1, freq = FALSE, 
     main=main,
     xlab=xlab, ylab = "Density", ...)
lines(density(qr1, adjust = 2),lwd=2)

qqnorm(qr1)
abline(0, 1, col = "red")

qr2 <- qnorm(p2)
hist(qr2, freq = FALSE, 
     main=main2,
     xlab=xlab2, ylab = "Density", ...)
lines(density(qr2, adjust = 2),lwd=2)

qqnorm(qr2)
abline(0, 1, col = "red")



if(test == TRUE){

#print(ks.test(qr1, "pnorm"))
print(shapiro.test(qr1))
#print(ks.test(qr2, "pnorm"))
print(shapiro.test(qr2))

}


}


}

